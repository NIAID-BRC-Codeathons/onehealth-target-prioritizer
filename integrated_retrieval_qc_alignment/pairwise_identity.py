
from pathlib import Path
import csv

ALIGNMENT_DIR = Path("alignments")
OUTPUT_DIR = Path("identity_results")
OUTPUT_DIR.mkdir(exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "all_pairwise_identity_ascending.csv"


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


def percent_identity(seq1, seq2):
    """
    Calculate pairwise amino-acid identity from an existing MSA.

    Columns containing a gap in either sequence are excluded.
    Identity = identical residues / compared non-gap positions * 100
    """

    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences have different lengths.")

    identical = 0
    compared = 0

    for aa1, aa2 in zip(seq1, seq2):

        # Exclude positions containing gaps
        if aa1 == "-" or aa2 == "-":
            continue

        compared += 1

        if aa1.upper() == aa2.upper():
            identical += 1

    if compared == 0:
        return 0.0, identical, compared

    identity = (identical / compared) * 100

    return identity, identical, compared


all_results = []

alignment_files = sorted(
    ALIGNMENT_DIR.glob("*_aligned.fasta")
)

for alignment_file in alignment_files:

    protein = alignment_file.stem.replace("_aligned", "")

    records = read_fasta(alignment_file)

    print(
        f"Processing {protein}: "
        f"{len(records)} sequences"
    )

    # Every unique pair
    for i in range(len(records)):

        header1, seq1 = records[i]

        for j in range(i + 1, len(records)):

            header2, seq2 = records[j]

            identity, identical, compared = percent_identity(
                seq1,
                seq2
            )

            all_results.append({
                "protein": protein,
                "sequence_1": header1,
                "sequence_2": header2,
                "identical_residues": identical,
                "positions_compared": compared,
                "percent_identity": round(identity, 2),
            })


# Sort LOWEST → HIGHEST identity
all_results.sort(
    key=lambda row: row["percent_identity"]
)


with open(OUTPUT_FILE, "w", newline="") as handle:

    fields = [
        "protein",
        "sequence_1",
        "sequence_2",
        "identical_residues",
        "positions_compared",
        "percent_identity",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(all_results)


print()
print("=" * 70)
print("PAIRWISE IDENTITY ANALYSIS COMPLETE")
print("=" * 70)

print(
    f"Total pairwise comparisons: "
    f"{len(all_results)}"
)

if all_results:

    print()
    print("Lowest identity:")
    print(
        f"{all_results[0]['protein']} | "
        f"{all_results[0]['percent_identity']}%"
    )

    print()
    print("Highest identity:")
    print(
        f"{all_results[-1]['protein']} | "
        f"{all_results[-1]['percent_identity']}%"
    )

print()
print(f"CSV written to: {OUTPUT_FILE}")

