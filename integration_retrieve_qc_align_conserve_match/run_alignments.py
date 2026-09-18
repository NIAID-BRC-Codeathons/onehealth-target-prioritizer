#!/usr/bin/env python3

from pathlib import Path
import re
import subprocess

INPUT = Path.cwd() / "converted.fasta"
OUTDIR = Path("alignments")
OUTDIR.mkdir(exist_ok=True)

# Protein groups and annotation terms expected in FASTA headers
TARGETS = {
    "nsp5_3CLpro": [
        r"\bnsp5\b",
        r"3C-like proteinase",
        r"3CLpro",
        r"main protease",
    ],
    "nsp12_RdRp": [
        r"\bnsp12\b",
        r"RNA-dependent RNA polymerase",
        r"\bRdRp\b",
    ],
    "nsp13_helicase": [
        r"\bnsp13\b",
        r"\bhelicase\b",
    ],
    "nsp14_ExoN": [
        r"\bnsp14\b",
        r"3['’]?-to-5['’]? exonuclease",
        r"\bexonuclease\b",
        r"\bExoN\b",
    ],
    "nsp15_NendoU": [
        r"\bnsp15\b",
        r"endoRNAse",
        r"endoribonuclease",
        r"\bNendoU\b",
    ],
    "nsp16_2O_MTase": [
        r"\bnsp16\b",
        r"2['’]?-O-ribose methyltransferase",
        r"2['’]?-O-methyltransferase",
    ],
}


def read_fasta(path):
    records = []
    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))

                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if header is not None:
            records.append((header, "".join(sequence)))

    return records


def matches(header, patterns):
    return any(
        re.search(pattern, header, flags=re.IGNORECASE)
        for pattern in patterns
    )


def write_fasta(records, path):
    with open(path, "w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")

            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i + 80] + "\n")


records = read_fasta(INPUT)

print(f"Total sequences in input: {len(records)}")
print()

for group, patterns in TARGETS.items():

    selected = [
        record
        for record in records
        if matches(record[0], patterns)
    ]

    print(f"{group}: {len(selected)} sequences")

    if len(selected) < 2:
        print("  Skipping: fewer than 2 sequences.")
        continue

    input_file = OUTDIR / f"{group}.fasta"
    output_file = OUTDIR / f"{group}_aligned.fasta"

    write_fasta(selected, input_file)

    print(f"  Running MAFFT...")

    with open(output_file, "w") as output_handle:
        subprocess.run(
            [
                "mafft",
                "--auto",
                "--thread",
                "-1",
                str(input_file),
            ],
            stdout=output_handle,
            check=True,
        )

    print(f"  Alignment written to {output_file}")
    print()

print("All available alignments completed.")


