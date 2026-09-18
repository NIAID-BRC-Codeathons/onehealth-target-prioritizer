#!/usr/bin/env python3
"""
integrated_viral_pipeline.py

Run the viral-protein workflow end-to-end:

    viral_protein_pipeline.py fetch
        -> json_to_fasta.py
        -> run_alignments.py
        -> pairwise_identity.py
        -> analyze_conservation.py
        -> conserved_regions.py (for every *_aligned.fasta)
        -> cov_pos_ab.tsv peptide search

The runner creates a run-specific working directory, so outputs from one
run do not overwrite outputs from another run.

Usage
-----
    python integrated_viral_pipeline.py --taxid-csv taxonomy_ids.csv

Optional:
    --limit 200
    --host-filter
    --work-dir results/my_run
    --python /path/to/python
    --fasta-width 60
    --conserved-min-identity 0.95
    --conserved-min-occupancy 0.50
    --conserved-min-length 10
    --conserved-max-length 15
    --conserved-weighting none
    --conserved-gap-votes

Notes
-----
The supplied run_alignments.py currently hard-codes:
    ~/Downloads/converted.fasta

This runner automatically places/copies the generated FASTA at that location
for the alignment stage, so run_alignments.py can be used without modification.

The conservation-window stage uses conserved_regions.py and is run once
for every alignment matching alignments/*_aligned.fasta.
The resulting peptide FASTA files are then searched exactly (case-insensitive)
against the `peptide` column of cov_pos_ab.tsv.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import csv
from datetime import datetime
from pathlib import Path


def log(message: str) -> None:
    print(f"[PIPELINE] {message}", flush=True)


def run_command(
    command: list[str],
    cwd: Path,
    log_file: Path,
) -> None:
    """Run a command, stream output, and stop on failure."""
    log(f"Running: {' '.join(command)}")

    with log_file.open("a", encoding="utf-8") as log_handle:
        log_handle.write(
            f"\n\n===== COMMAND: {' '.join(command)} =====\n"
        )

        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdout is not None

        for line in process.stdout:
            print(line, end="")
            log_handle.write(line)

        return_code = process.wait()

        log_handle.write(
            f"\n===== EXIT CODE: {return_code} =====\n"
        )

    if return_code != 0:
        raise RuntimeError(
            f"Command failed with exit code {return_code}: "
            f"{' '.join(command)}"
        )


def require_file(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"{description} was not created: {path}"
        )

    if not path.is_file():
        raise FileNotFoundError(
            f"{description} is not a file: {path}"
        )


def require_dir(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"{description} was not created: {path}"
        )

    if not path.is_dir():
        raise FileNotFoundError(
            f"{description} is not a directory: {path}"
        )


def count_json_objects(path: Path) -> int:
    """
    Count concatenated JSON objects using JSONDecoder.raw_decode().
    This matches the format produced by viral_protein_pipeline.py when
    multiple Taxonomy IDs are supplied.
    """
    decoder = json.JSONDecoder()
    text = path.read_text(encoding="utf-8")
    position = 0
    count = 0

    while position < len(text):
        while position < len(text) and text[position].isspace():
            position += 1

        if position >= len(text):
            break

        _, position = decoder.raw_decode(text, position)
        count += 1

    return count


def find_alignment_outputs(alignment_dir: Path) -> list[Path]:
    return sorted(alignment_dir.glob("*_aligned.fasta"))


def read_cov_peptides(path: Path) -> list[dict[str, str]]:
    """Load cov_pos_ab.tsv rows for substring peptide searching."""

    if not path.exists():
        raise FileNotFoundError(f"cov_pos_ab.tsv not found: {path}")

    rows_out: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline=""
    ) as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames or "peptide" not in reader.fieldnames:
            raise ValueError(
                f"{path} does not contain a 'peptide' column"
            )

        for row in reader:
            peptide = (row.get("peptide") or "").strip().upper()

            if not peptide:
                continue

            # Store the normalized peptide back into the row.
            row["peptide"] = peptide
            rows_out.append(row)

    return rows_out


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Read FASTA records as (header, sequence) pairs."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq).upper()))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq).upper()))

    return records


def read_matrix_rows(path: Path) -> dict[str, dict[str, str]]:
    """Read conserved-region matrix rows keyed by sequence ID."""
    if not path.exists():
        return {}

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if len(rows) < 4:
        return {}

    # The matrix has two header rows followed by a consensus row.
    field_row = rows[1]
    result: dict[str, dict[str, str]] = {}

    for row in rows[3:]:
        if not row or not row[0]:
            continue
        values: dict[str, str] = {
            "sequence_id": row[0],
        }
        # Preserve the available matrix columns without assuming one stretch.
        for i in range(1, min(len(row), len(field_row))):
            values[f"column_{i}"] = row[i]
        result[row[0]] = values

    return result


def search_cov_pos_ab(
    conserved_regions_dir: Path,
    cov_pos_ab_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Search every conserved-region peptide against cov_pos_ab.tsv."""
    fasta_files = sorted(conserved_regions_dir.glob("*_conserved_regions.fasta"))
    if not fasta_files:
        raise RuntimeError(
            f"No *_conserved_regions.fasta files found in {conserved_regions_dir}"
        )

    peptide_lookup = read_cov_peptides(cov_pos_ab_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    total_candidates = 0
    total_matches = 0

    for fasta_file in fasta_files:
        stem = fasta_file.name.replace("_conserved_regions.fasta", "")
        matrix_file = conserved_regions_dir / f"{stem}_matrix.tsv"
        matrix_rows = read_matrix_rows(matrix_file)

        out_file = output_dir / f"{stem}_cov_pos_ab_matches.tsv"
        fieldnames = [
            "source_fasta",
            "sequence_id",
            "peptide",
            "peptide_length",
            "region",
            "residues",
            "identity",
            "matched",
            "structure_id",
            "source_organism",
            "protein",
            "protein_all",
            "receptor",
            "receptor_confidence",
            "is_bcell",
            "is_tcell",
            "protein_role",
            "protein_receptor",
        ]

        rows_out: list[dict[str, object]] = []

        for header, sequence in read_fasta_records(fasta_file):
            sequence_id = header.split()[0]
            peptide = sequence.upper()
            total_candidates += 1

            matches = [
                row
                for row in peptide_lookup
                if peptide in row.get("peptide", "")
            ]

            matrix_row = matrix_rows.get(sequence_id, {})
            # Header contains region/residues/identity fields for the generated FASTA.
            header_fields: dict[str, str] = {}
            for token in header.split():
                if "=" in token:
                    key, value = token.split("=", 1)
                    header_fields[key] = value

            if matches:
                total_matches += 1
                for match in matches:
                    rows_out.append({
                        "source_fasta": fasta_file.name,
                        "sequence_id": sequence_id,
                        "peptide": peptide,
                        "peptide_length": len(peptide),
                        "region": header_fields.get("region", ""),
                        "residues": header_fields.get("residues", ""),
                        "identity": header_fields.get("identity", ""),
                        "matched": "1",
                        "structure_id": match.get("structure_id", ""),
                        "source_organism": match.get("source_organism", ""),
                        "protein": match.get("protein", ""),
                        "protein_all": match.get("protein_all", ""),
                        "receptor": match.get("receptor", ""),
                        "receptor_confidence": match.get("receptor_confidence", ""),
                        "is_bcell": match.get("is_bcell", ""),
                        "is_tcell": match.get("is_tcell", ""),
                        "protein_role": match.get("protein_role", ""),
                        "protein_receptor": match.get("protein_receptor", ""),
                    })
            else:
                rows_out.append({
                    "source_fasta": fasta_file.name,
                    "sequence_id": sequence_id,
                    "peptide": peptide,
                    "peptide_length": len(peptide),
                    "region": header_fields.get("region", ""),
                    "residues": header_fields.get("residues", ""),
                    "identity": header_fields.get("identity", ""),
                    "matched": "0",
                    "structure_id": "",
                    "source_organism": "",
                    "protein": "",
                    "protein_all": "",
                    "receptor": "",
                    "receptor_confidence": "",
                    "is_bcell": "",
                    "is_tcell": "",
                    "protein_role": "",
                    "protein_receptor": "",
                })

        with out_file.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows_out)

        summary_rows.append({
            "source_fasta": fasta_file.name,
            "candidates": len(read_fasta_records(fasta_file)),
            "matched_candidates": sum(1 for r in rows_out if r["matched"] == "1"),
            "output": str(out_file),
        })

    summary_file = output_dir / "cov_pos_ab_match_summary.tsv"
    with summary_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_fasta", "candidates", "matched_candidates", "output"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    return {
        "cov_pos_ab": str(cov_pos_ab_path),
        "source_files": len(fasta_files),
        "candidate_peptides": total_candidates,
        "matched_candidate_peptides": total_matches,
        "summary": str(summary_file),
        "files": summary_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run fetch -> JSON-to-FASTA -> MAFFT alignments -> "
            "pairwise identity -> conservation -> cov_pos_ab peptide search."
        )
    )

    parser.add_argument(
        "--taxid-csv",
        required=True,
        help=(
            "CSV file containing a column named 'Taxonomy ID'. "
            "Passed to viral_protein_pipeline.py fetch."
        ),
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=200,
        help="Maximum number of records requested per Taxonomy ID (default: 200).",
    )

    parser.add_argument(
        "--host-filter",
        action="store_true",
        help="Enable the host filter in viral_protein_pipeline.py.",
    )

    parser.add_argument(
        "--work-dir",
        default=None,
        help=(
            "Directory for this pipeline run. "
            "Default: ./pipeline_runs/YYYYMMDD_HHMMSS"
        ),
    )

    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used for all stages (default: current interpreter).",
    )

    parser.add_argument(
        "--fasta-width",
        type=int,
        default=60,
        help="FASTA line width passed to json_to_fasta.py (default: 60).",
    )

    parser.add_argument(
        "--cov-pos-ab",
        default="cov_pos_ab.tsv",
        help=(
            "Path to cov_pos_ab.tsv. The file must contain a 'peptide' column. "
            "Default: cov_pos_ab.tsv in the project directory."
        ),
    )

    parser.add_argument(
        "--keep-existing-run-dir",
        action="store_true",
        help=(
            "Allow an existing work directory to be reused. "
            "Without this option, an existing directory causes an error."
        ),
    )

    parser.add_argument(
        "--conserved-min-identity",
        type=float,
        default=0.95,
        help="Minimum conserved-column identity; accepts 0.95 or 95 (default: 0.95).",
    )

    parser.add_argument(
        "--conserved-min-occupancy",
        type=float,
        default=0.50,
        help="Minimum non-gap occupancy for a column (default: 0.50).",
    )

    parser.add_argument(
        "--conserved-min-length",
        type=int,
        default=10,
        help="Shortest peptide window to report (default: 10).",
    )

    parser.add_argument(
        "--conserved-max-length",
        type=int,
        default=15,
        help="Longest peptide window to report (default: 15).",
    )

    parser.add_argument(
        "--conserved-weighting",
        choices=["none", "henikoff"],
        default="none",
        help="Sequence weighting for conserved-region analysis (default: none).",
    )

    parser.add_argument(
        "--conserved-gap-votes",
        action="store_true",
        help="Count gaps in the conserved-column identity denominator.",
    )

    args = parser.parse_args()

    taxid_csv = Path(args.taxid_csv).resolve()
    cov_pos_ab_path = Path(args.cov_pos_ab)

    if not taxid_csv.exists():
        print(f"ERROR: Taxonomy CSV not found: {taxid_csv}", file=sys.stderr)
        return 2

    if args.limit <= 0:
        print("ERROR: --limit must be greater than 0.", file=sys.stderr)
        return 2

    if args.fasta_width <= 0:
        print("ERROR: --fasta-width must be greater than 0.", file=sys.stderr)
        return 2

    def normalize_fraction(value: float, name: str) -> float:
        # Accept either 0.95 or 95.
        normalized = value / 100.0 if value > 1.0 else value
        if not 0.0 <= normalized <= 1.0:
            print(
                f"ERROR: {name} must be between 0 and 1, or 0 and 100 as a percentage.",
                file=sys.stderr,
            )
            raise ValueError(name)
        return normalized

    try:
        conserved_min_identity = normalize_fraction(
            args.conserved_min_identity, "--conserved-min-identity"
        )
        conserved_min_occupancy = normalize_fraction(
            args.conserved_min_occupancy, "--conserved-min-occupancy"
        )
    except ValueError:
        return 2

    if args.conserved_min_length < 1:
        print("ERROR: --conserved-min-length must be >= 1.", file=sys.stderr)
        return 2

    if args.conserved_max_length < args.conserved_min_length:
        print(
            "ERROR: --conserved-max-length must be >= --conserved-min-length.",
            file=sys.stderr,
        )
        return 2

    project_dir = Path(__file__).resolve().parent
    if not cov_pos_ab_path.is_absolute():
        cov_pos_ab_path = (project_dir / cov_pos_ab_path).resolve()

    viral_script = project_dir / "viral_protein_pipeline.py"
    fasta_script = project_dir / "json_to_fasta.py"
    alignment_script = project_dir / "run_alignments.py"
    identity_script = project_dir / "pairwise_identity.py"
    conservation_script = project_dir / "analyze_conservation.py"
    conserved_regions_script = project_dir / "conserved_regions.py"

    required_scripts = [
        viral_script,
        fasta_script,
        alignment_script,
        identity_script,
        conservation_script,
        conserved_regions_script,
    ]

    for script in required_scripts:
        if not script.exists():
            print(
                f"ERROR: Required script not found: {script}",
                file=sys.stderr,
            )
            return 2

    if args.work_dir:
        work_dir = Path(args.work_dir)
        if not work_dir.is_absolute():
            work_dir = project_dir / work_dir
        work_dir = work_dir.resolve()
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        work_dir = project_dir / "pipeline_runs" / timestamp

    if work_dir.exists() and not args.keep_existing_run_dir:
        print(
            f"ERROR: Work directory already exists: {work_dir}\n"
            "Use --keep-existing-run-dir to reuse it, or choose another directory.",
            file=sys.stderr,
        )
        return 2

    work_dir.mkdir(parents=True, exist_ok=True)

    # Stage-specific paths
    json_output = work_dir / "retrieved_proteins.json"
    fasta_output = work_dir / "converted.fasta"
    alignment_dir = work_dir / "alignments"
    identity_dir = work_dir / "identity_results"
    conservation_dir = work_dir / "conservation_results"
    conserved_regions_dir = work_dir / "conserved_regions"
    log_file = work_dir / "pipeline.log"

    # run_alignments.py uses a hard-coded path under ~/Downloads.
    # downloads_fasta = Path.home() / "Downloads" / "converted.fasta"

    try:
        log("Starting integrated viral protein pipeline.")
        log(f"Work directory: {work_dir}")
        log(f"Taxonomy CSV: {taxid_csv}")

        # -----------------------------------------------------------------
        # Stage 1: Fetch
        # -----------------------------------------------------------------
        fetch_command = [
            args.python,
            str(viral_script),
            "fetch",
            "--taxid",
            str(taxid_csv),
            "--limit",
            str(args.limit),
            "--output",
            str(json_output),
        ]

        if args.host_filter:
            fetch_command.append("--host-filter")

        run_command(fetch_command, cwd=work_dir, log_file=log_file)
        require_file(json_output, "Fetch JSON output")

        object_count = count_json_objects(json_output)
        log(
            f"Fetch stage complete: {object_count} JSON object(s) "
            f"written to {json_output}"
        )

        # -----------------------------------------------------------------
        # Stage 2: JSON -> FASTA
        # -----------------------------------------------------------------
        fasta_command = [
            args.python,
            str(fasta_script),
            str(json_output),
            "--output",
            str(fasta_output),
            "--width",
            str(args.fasta_width),
        ]

        run_command(fasta_command, cwd=work_dir, log_file=log_file)
        require_file(fasta_output, "FASTA output")

        # -----------------------------------------------------------------
        # Stage 3: Prepare hard-coded path expected by run_alignments.py
        # -----------------------------------------------------------------
        # downloads_fasta.parent.mkdir(parents=True, exist_ok=True)

        # Back up an existing converted.fasta rather than silently deleting it.
        # backup_path = None
        # if downloads_fasta.exists():
        #     timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        #     backup_path = downloads_fasta.with_name(
        #         f"converted.fasta.pipeline_backup_{timestamp}"
        #     )
        #     shutil.copy2(downloads_fasta, backup_path)
        #     log(f"Backed up existing {downloads_fasta} to {backup_path}")

        # shutil.copy2(fasta_output, downloads_fasta)
        # log(
        #     "Copied generated FASTA to the path expected by run_alignments.py: "
        #     f"{downloads_fasta}"
        # )

        # -----------------------------------------------------------------
        # Stage 4: Run MAFFT alignments
        # -----------------------------------------------------------------
        # run_alignments.py creates "alignments" relative to its cwd.
        run_command(
            [args.python, str(alignment_script)],
            cwd=work_dir,
            log_file=log_file,
        )

        require_dir(alignment_dir, "Alignment directory")

        alignment_files = find_alignment_outputs(alignment_dir)

        if not alignment_files:
            raise RuntimeError(
                "Alignment stage completed but no *_aligned.fasta files were produced."
            )

        log(
            f"Alignment stage complete: {len(alignment_files)} aligned FASTA file(s)."
        )

        # -----------------------------------------------------------------
        # Stage 5: Pairwise identity
        # -----------------------------------------------------------------
        run_command(
            [args.python, str(identity_script)],
            cwd=work_dir,
            log_file=log_file,
        )

        require_file(
            identity_dir / "all_pairwise_identity_ascending.csv",
            "Pairwise identity output",
        )

        log("Pairwise identity stage complete.")

        # -----------------------------------------------------------------
        # Stage 6: Conservation analysis
        # -----------------------------------------------------------------
        run_command(
            [args.python, str(conservation_script)],
            cwd=work_dir,
            log_file=log_file,
        )

        require_file(
            conservation_dir / "conservation_summary.csv",
            "Conservation summary",
        )

        log("Conservation analysis stage complete.")

        # -----------------------------------------------------------------
        # Stage 7: Conserved peptide-region extraction
        # -----------------------------------------------------------------
        conserved_regions_dir.mkdir(parents=True, exist_ok=True)

        conserved_outputs = []
        for alignment_file in alignment_files:
            stem = alignment_file.stem  # e.g. nsp5_3CLpro_aligned
            matrix_output = conserved_regions_dir / f"{stem}_matrix.tsv"
            fasta_output_regions = conserved_regions_dir / f"{stem}_conserved_regions.fasta"

            conserved_command = [
                args.python,
                str(conserved_regions_script),
                str(alignment_file),
                "--min-identity",
                str(conserved_min_identity),
                "--min-occupancy",
                str(conserved_min_occupancy),
                "--min-length",
                str(args.conserved_min_length),
                "--max-length",
                str(args.conserved_max_length),
                "--weighting",
                args.conserved_weighting,
                "--matrix-out",
                str(matrix_output),
                "--fasta-out",
                str(fasta_output_regions),
            ]

            if args.conserved_gap_votes:
                conserved_command.append("--gap-votes")

            run_command(
                conserved_command,
                cwd=work_dir,
                log_file=log_file,
            )

            require_file(matrix_output, f"Conserved-region matrix for {alignment_file.name}")
            require_file(fasta_output_regions, f"Conserved-region FASTA for {alignment_file.name}")

            conserved_outputs.append({
                "alignment": str(alignment_file),
                "matrix": str(matrix_output),
                "fasta": str(fasta_output_regions),
            })

            log(
                f"Conserved-region extraction complete for {alignment_file.name}: "
                f"{matrix_output.name}, {fasta_output_regions.name}"
            )

        if not conserved_outputs:
            raise RuntimeError("No conserved-region outputs were produced.")

        log(
            f"Conserved-region stage complete: {len(conserved_outputs)} alignment file(s) analyzed."
        )

        # -----------------------------------------------------------------
        # Stage 8: Search conserved peptides in cov_pos_ab.tsv
        # -----------------------------------------------------------------
        if not cov_pos_ab_path.exists():
            raise FileNotFoundError(f"cov_pos_ab.tsv not found: {cov_pos_ab_path}")

        cov_match_dir = work_dir / "cov_pos_ab_matches"
        cov_match_results = search_cov_pos_ab(
            conserved_regions_dir=conserved_regions_dir,
            cov_pos_ab_path=cov_pos_ab_path,
            output_dir=cov_match_dir,
        )

        log(
            "cov_pos_ab peptide search complete: "
            f"{cov_match_results['matched_candidate_peptides']} / "
            f"{cov_match_results['candidate_peptides']} candidate peptides matched."
        )

        # -----------------------------------------------------------------
        # Manifest
        # -----------------------------------------------------------------
        manifest = {
            "taxid_csv": str(taxid_csv),
            "limit": args.limit,
            "host_filter": args.host_filter,
            "fasta_width": args.fasta_width,
            "conserved_region_settings": {
                "min_identity": conserved_min_identity,
                "min_occupancy": conserved_min_occupancy,
                "min_length": args.conserved_min_length,
                "max_length": args.conserved_max_length,
                "weighting": args.conserved_weighting,
                "gap_votes": args.conserved_gap_votes,
            },
            "work_dir": str(work_dir),
            "outputs": {
                "json": str(json_output),
                "fasta": str(fasta_output),
                "alignments": [str(p) for p in alignment_files],
                "pairwise_identity": str(
                    identity_dir / "all_pairwise_identity_ascending.csv"
                ),
                "conservation_summary": str(
                    conservation_dir / "conservation_summary.csv"
                ),
                "conservation_dir": str(conservation_dir),
                "conserved_regions": conserved_outputs,
                "conserved_regions_dir": str(conserved_regions_dir),
                "cov_pos_ab": cov_match_results,
                "cov_pos_ab_matches_dir": str(cov_match_dir),
            },
        }

        manifest_path = work_dir / "pipeline_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2),
            encoding="utf-8",
        )

        log("Pipeline completed successfully.")
        log(f"Results: {work_dir}")
        log(f"Manifest: {manifest_path}")

        # Restore pre-existing ~/Downloads/converted.fasta when applicable.
        # if backup_path is not None:
        #     shutil.copy2(backup_path, downloads_fasta)
        #     backup_path.unlink(missing_ok=True)
        #     log(f"Restored original {downloads_fasta}")
        # else:
        #     # Do not remove a file that did not exist before the pipeline.
        #     downloads_fasta.unlink(missing_ok=True)
        #     log(f"Removed temporary {downloads_fasta}")

        # return 0

    except Exception as exc:
        log(f"PIPELINE FAILED: {exc}")
        log(f"Pipeline log: {log_file}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
