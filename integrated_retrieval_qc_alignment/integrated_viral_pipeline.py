#!/usr/bin/env python3
"""
integrated_viral_pipeline.py

Run the viral-protein workflow end-to-end:

    viral_protein_pipeline.py fetch
        -> json_to_fasta.py
        -> run_alignments.py
        -> pairwise_identity.py
        -> analyze_conservation.py

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

Notes
-----
The supplied run_alignments.py currently hard-codes:
    ~/Downloads/converted.fasta

This runner automatically places/copies the generated FASTA at that location
for the alignment stage, so run_alignments.py can be used without modification.

The final script supplied by the user is named analyze_conservation.py.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run fetch -> JSON-to-FASTA -> MAFFT alignments -> "
            "pairwise identity -> conservation analysis."
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
        "--keep-existing-run-dir",
        action="store_true",
        help=(
            "Allow an existing work directory to be reused. "
            "Without this option, an existing directory causes an error."
        ),
    )

    args = parser.parse_args()

    taxid_csv = Path(args.taxid_csv).resolve()

    if not taxid_csv.exists():
        print(f"ERROR: Taxonomy CSV not found: {taxid_csv}", file=sys.stderr)
        return 2

    if args.limit <= 0:
        print("ERROR: --limit must be greater than 0.", file=sys.stderr)
        return 2

    if args.fasta_width <= 0:
        print("ERROR: --fasta-width must be greater than 0.", file=sys.stderr)
        return 2

    project_dir = Path(__file__).resolve().parent

    viral_script = project_dir / "viral_protein_pipeline.py"
    fasta_script = project_dir / "json_to_fasta.py"
    alignment_script = project_dir / "run_alignments.py"
    identity_script = project_dir / "pairwise_identity.py"
    conservation_script = project_dir / "analyze_conservation.py"

    required_scripts = [
        viral_script,
        fasta_script,
        alignment_script,
        identity_script,
        conservation_script,
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
    log_file = work_dir / "pipeline.log"

    # run_alignments.py uses a hard-coded path under ~/Downloads.
    downloads_fasta = Path.home() / "Downloads" / "converted.fasta"

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
        downloads_fasta.parent.mkdir(parents=True, exist_ok=True)

        # Back up an existing converted.fasta rather than silently deleting it.
        backup_path = None
        if downloads_fasta.exists():
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_path = downloads_fasta.with_name(
                f"converted.fasta.pipeline_backup_{timestamp}"
            )
            shutil.copy2(downloads_fasta, backup_path)
            log(f"Backed up existing {downloads_fasta} to {backup_path}")

        shutil.copy2(fasta_output, downloads_fasta)
        log(
            "Copied generated FASTA to the path expected by run_alignments.py: "
            f"{downloads_fasta}"
        )

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
        # Manifest
        # -----------------------------------------------------------------
        manifest = {
            "taxid_csv": str(taxid_csv),
            "limit": args.limit,
            "host_filter": args.host_filter,
            "fasta_width": args.fasta_width,
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
        if backup_path is not None:
            shutil.copy2(backup_path, downloads_fasta)
            backup_path.unlink(missing_ok=True)
            log(f"Restored original {downloads_fasta}")
        else:
            # Do not remove a file that did not exist before the pipeline.
            downloads_fasta.unlink(missing_ok=True)
            log(f"Removed temporary {downloads_fasta}")

        return 0

    except Exception as exc:
        log(f"PIPELINE FAILED: {exc}")
        log(f"Pipeline log: {log_file}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
