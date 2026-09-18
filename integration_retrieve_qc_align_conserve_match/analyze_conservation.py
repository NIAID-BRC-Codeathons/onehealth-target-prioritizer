#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv

ALIGNMENT_DIR = Path("alignments")
OUTPUT_DIR = Path("conservation_results")
OUTPUT_DIR.mkdir(exist_ok=True)

# Threshold for calling a site highly conserved
CONSERVATION_THRESHOLD = 0.90


def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq)))

    return records


def analyze_alignment(path):

    records = read_fasta(path)

    if not records:
        return None

    lengths = {len(seq) for _, seq in records}

    if len(lengths) != 1:
        raise ValueError(
            f"{path}: sequences do not have equal aligned lengths."
        )

    alignment_length = lengths.pop()
    n_sequences = len(records)

    results = []

    for pos in range(alignment_length):

        column = [seq[pos] for _, seq in records]

        gaps = column.count("-")
        residues = [aa for aa in column if aa != "-"]

        if residues:
            counts = Counter(residues)
            consensus, consensus_count = counts.most_common(1)[0]

            # Conservation is calculated against all sequences.
            conservation = consensus_count / n_sequences

        else:
            consensus = "-"
            conservation = 0

        gap_fraction = gaps / n_sequences

        identical = (
            len(set(column)) == 1
            and "-" not in column
        )

        results.append({
            "alignment_position": pos + 1,
            "consensus": consensus,
            "conservation": conservation,
            "gap_fraction": gap_fraction,
            "identical": identical,
        })

    return records, results


def conserved_blocks(results, threshold=0.90):

    blocks = []
    start = None

    for r in results:

        conserved = (
            r["conservation"] >= threshold
            and r["gap_fraction"] == 0
        )

        if conserved and start is None:
            start = r["alignment_position"]

        elif not conserved and start is not None:
            end = r["alignment_position"] - 1

            if end - start + 1 >= 3:
                blocks.append((start, end))

            start = None

    if start is not None:
        end = results[-1]["alignment_position"]

        if end - start + 1 >= 3:
            blocks.append((start, end))

    return blocks


summary = []

files = sorted(
    ALIGNMENT_DIR.glob("*_aligned.fasta")
)

for path in files:

    protein = path.stem.replace("_aligned", "")

    analyzed = analyze_alignment(path)

    if analyzed is None:
        continue

    records, results = analyzed

    n_sequences = len(records)
    alignment_length = len(results)

    identical_sites = sum(
        r["identical"] for r in results
    )

    conserved_sites = sum(
        r["conservation"] >= CONSERVATION_THRESHOLD
        and r["gap_fraction"] == 0
        for r in results
    )

    blocks = conserved_blocks(
        results,
        CONSERVATION_THRESHOLD
    )

    # Write per-position conservation table
    csv_file = OUTPUT_DIR / f"{protein}_conservation.csv"

    with open(csv_file, "w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "alignment_position",
                "consensus",
                "conservation",
                "gap_fraction",
                "identical",
            ],
        )

        writer.writeheader()
        writer.writerows(results)

    # Write conserved blocks
    block_file = OUTPUT_DIR / f"{protein}_conserved_blocks.txt"

    with open(block_file, "w") as handle:

        for start, end in blocks:

            sequence = "".join(
                results[i - 1]["consensus"]
                for i in range(start, end + 1)
            )

            handle.write(
                f"{start}-{end}\t"
                f"length={end-start+1}\t"
                f"{sequence}\n"
            )

    summary.append({
        "protein": protein,
        "sequences": n_sequences,
        "alignment_length": alignment_length,
        "identical_sites": identical_sites,
        "highly_conserved_sites": conserved_sites,
        "percent_identical_sites":
            round(100 * identical_sites / alignment_length, 2),
        "percent_highly_conserved":
            round(100 * conserved_sites / alignment_length, 2),
        "conserved_blocks": len(blocks),
    })


summary_file = OUTPUT_DIR / "conservation_summary.csv"

with open(summary_file, "w", newline="") as handle:

    fields = [
        "protein",
        "sequences",
        "alignment_length",
        "identical_sites",
        "highly_conserved_sites",
        "percent_identical_sites",
        "percent_highly_conserved",
        "conserved_blocks",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(summary)


print("\nCONSERVATION SUMMARY")
print("=" * 95)

for row in summary:

    print(
        f"{row['protein']:20s} "
        f"N={row['sequences']:3d}  "
        f"Length={row['alignment_length']:4d}  "
        f"Identical={row['identical_sites']:4d} "
        f"({row['percent_identical_sites']:6.2f}%)  "
        f">=90%={row['highly_conserved_sites']:4d} "
        f"({row['percent_highly_conserved']:6.2f}%)  "
        f"Blocks={row['conserved_blocks']}"
    )

print()
print(f"Results written to: {OUTPUT_DIR}/")
