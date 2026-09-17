#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


def clean_sequence(sequence):
    """Remove whitespace and convert sequence to uppercase."""
    return re.sub(r"\s+", "", str(sequence)).upper()


def wrap_sequence(sequence, width=60):
    """Wrap a sequence to fixed-width lines."""
    for i in range(0, len(sequence), width):
        yield sequence[i:i + width]


def fasta_header(record):
    """Create a FASTA header."""
    accession = record.get("accession", "unknown")
    title = record.get("title", "unknown protein")
    organism = record.get("organism", "")

    if organism:
        return f">{accession} {title} | {organism}"
    return f">{accession} {title}"


def read_concatenated_json(filename):
    """
    Read a file containing multiple JSON objects concatenated together.

    Example:

        {...}
        {...}
        {...}

    Returns each JSON object one at a time.
    """
    decoder = json.JSONDecoder()

    with open(filename, "r", encoding="utf-8") as f:
        text = f.read()

    position = 0
    length = len(text)

    while position < length:
        # Skip whitespace between JSON objects
        while position < length and text[position].isspace():
            position += 1

        if position >= length:
            break

        try:
            obj, next_position = decoder.raw_decode(text, position)
        except json.JSONDecodeError as e:
            raise ValueError(
                f"Could not parse JSON near character {e.pos} "
                f"(line {e.lineno}, column {e.colno}): {e.msg}"
            ) from e

        yield obj
        position = next_position


def convert_to_fasta(input_file, output_file, width=60):
    """Convert all sequence records in concatenated JSON to FASTA."""

    total_objects = 0
    total_records = 0
    skipped_records = 0
    length_warnings = 0

    with open(output_file, "w", encoding="utf-8") as out:

        for json_object in read_concatenated_json(input_file):
            total_objects += 1

            if not isinstance(json_object, dict):
                print(
                    f"Skipping JSON object {total_objects}: "
                    f"expected an object, got {type(json_object).__name__}"
                )
                continue

            records = json_object.get("records", [])

            if not isinstance(records, list):
                print(
                    f"Warning: object {total_objects} has a non-list "
                    f"'records' field; skipping it."
                )
                continue

            for record in records:

                if not isinstance(record, dict):
                    skipped_records += 1
                    continue

                accession = record.get("accession")
                sequence = record.get("sequence")

                if not accession or not sequence:
                    print(
                        "Skipping record missing accession or sequence:",
                        record
                    )
                    skipped_records += 1
                    continue

                sequence = clean_sequence(sequence)

                # Check reported length against actual sequence length
                reported_length = record.get("length")

                if reported_length is not None:
                    try:
                        reported_length = int(reported_length)

                        if len(sequence) != reported_length:
                            print(
                                f"WARNING: {accession}: "
                                f"reported length={reported_length}, "
                                f"actual length={len(sequence)}"
                            )
                            length_warnings += 1

                    except (ValueError, TypeError):
                        pass

                # Write FASTA entry
                out.write(fasta_header(record) + "\n")

                for line in wrap_sequence(sequence, width):
                    out.write(line + "\n")

                total_records += 1

    print("\nConversion complete")
    print("-------------------")
    print(f"JSON objects read : {total_objects}")
    print(f"FASTA records     : {total_records}")
    print(f"Records skipped   : {skipped_records}")
    print(f"Length warnings   : {length_warnings}")
    print(f"Output            : {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert concatenated JSON protein records to FASTA."
    )

    parser.add_argument(
        "input",
        help="Input JSON/text file containing concatenated JSON objects"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="converted.fasta",
        help="Output FASTA file (default: converted.fasta)"
    )

    parser.add_argument(
        "-w",
        "--width",
        type=int,
        default=60,
        help="Residues per FASTA line (default: 60)"
    )

    args = parser.parse_args()

    input_file = Path(args.input)

    if not input_file.exists():
        raise SystemExit(f"Input file not found: {input_file}")

    if args.width <= 0:
        raise SystemExit("--width must be greater than zero")

    try:
        convert_to_fasta(
            input_file=input_file,
            output_file=args.output,
            width=args.width
        )
    except Exception as e:
        raise SystemExit(f"Error: {e}")


if __name__ == "__main__":
    main()