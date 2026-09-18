# Integrated Viral Protein Analysis Pipeline

`integrated_viral_pipeline.py` runs an end-to-end viral protein workflow from NCBI protein retrieval through alignment, conservation analysis, conserved-peptide extraction, and exact peptide lookup against `cov_pos_ab.tsv`.

## Workflow

```text
Taxonomy ID CSV
      |
      v
viral_protein_pipeline.py fetch
      |
      v
retrieved_proteins.json
      |
      v
json_to_fasta.py
      |
      v
converted.fasta
      |
      v
run_alignments.py + MAFFT
      |
      v
alignments/*_aligned.fasta
      |
      +----------------------+-------------------------+
      |                      |                         |
      v                      v                         v
pairwise_identity.py  analyze_conservation.py  conserved_regions.py
      |                      |                         |
      v                      v                         +--> *_matrix.tsv
identity_results/     conservation_results/            +--> *_conserved_regions.fasta
                                                        |
                                                        v
                                                exact peptide lookup
                                                in cov_pos_ab.tsv
                                                        |
                                                        v
                                                cov_pos_ab_matches/
```

The runner uses a run-specific working directory and stops when a required stage fails or an expected output is missing.

## Required project files

Keep these files in the same directory as `integrated_viral_pipeline.py`:

```text
integrated_viral_pipeline.py
viral_protein_pipeline.py
json_to_fasta.py
run_alignments.py
pairwise_identity.py
analyze_conservation.py
conserved_regions.py
cov_pos_ab.tsv
taxonomy_ids.csv
```

`cov_pos_ab.tsv` may be stored elsewhere if its path is supplied with `--cov-pos-ab`.

## Dependencies

### Python

Python 3.9 or newer is recommended.

Check the installed version:

```bash
python --version
```

### Python packages

The workflow requires:

```text
pandas
numpy
biopython
```

Install them with:

```bash
python -m pip install pandas numpy biopython
```

A minimal `requirements.txt` is:

```text
pandas>=2.0
numpy
biopython
```

`json_to_fasta.py`, `pairwise_identity.py`, `analyze_conservation.py`, and the integrated runner otherwise rely primarily on Python standard-library modules.

### MAFFT

MAFFT is an external command-line dependency used by `run_alignments.py`. It must be installed separately and available on `PATH`.

Verify:

```bash
mafft --version
```

macOS/Homebrew:

```bash
brew install mafft
```

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install mafft
```

Conda/Bioconda:

```bash
conda install -c bioconda mafft
```

### Internet access

The fetch stage requires internet access to NCBI E-utilities.

## Taxonomy ID input

The required input CSV must contain a column named exactly:

```text
Taxonomy ID
```

Example `taxonomy_ids.csv`:

```csv
Taxonomy ID
2697049
694009
11137
```

Multiple Taxonomy IDs are supported.

## Basic usage

With `cov_pos_ab.tsv` in the same directory as the integrated runner:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv
```

To explicitly provide the peptide database:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --cov-pos-ab /path/to/cov_pos_ab.tsv
```

Example with additional settings:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --limit 500 \
    --host-filter \
    --work-dir results/run1 \
    --cov-pos-ab cov_pos_ab.tsv \
    --conserved-min-identity 0.95 \
    --conserved-min-occupancy 0.50 \
    --conserved-min-length 10 \
    --conserved-max-length 15
```

## Main command-line arguments

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--taxid-csv` | Yes | — | CSV containing the `Taxonomy ID` column. |
| `--limit` | No | `200` | Maximum NCBI records requested per Taxonomy ID. |
| `--host-filter` | No | Off | Enables the host filter in the retrieval script. |
| `--work-dir` | No | `pipeline_runs/YYYYMMDD_HHMMSS` | Output directory for the run. |
| `--python` | No | Current interpreter | Python executable used for component scripts. |
| `--fasta-width` | No | `60` | FASTA line width used by `json_to_fasta.py`. |
| `--cov-pos-ab` | No | `cov_pos_ab.tsv` | Peptide database TSV. Relative paths are resolved from the project directory. |
| `--keep-existing-run-dir` | No | Off | Allows reuse of an existing work directory. |

### Conserved-region arguments

| Argument | Default | Description |
| --- | ---: | --- |
| `--conserved-min-identity` | `0.95` | Minimum conserved-column identity. Values such as `95` are also accepted. |
| `--conserved-min-occupancy` | `0.50` | Minimum non-gap column occupancy. |
| `--conserved-min-length` | `10` | Shortest peptide window. |
| `--conserved-max-length` | `15` | Longest peptide window. |
| `--conserved-weighting` | `none` | `none` or `henikoff`. |
| `--conserved-gap-votes` | Off | Counts gaps in the identity denominator when enabled. |

Display all options with:

```bash
python integrated_viral_pipeline.py --help
```

## Stage 1 — Retrieve NCBI protein records

The runner executes `viral_protein_pipeline.py fetch`, passing the Taxonomy ID CSV, record limit, optional host filter, and JSON output path.

The retrieval stage stores metadata including accession, title, organism, length, and protein sequence in:

```text
retrieved_proteins.json
```

When multiple Taxonomy IDs are supplied, the retrieval function appends one JSON object per Taxonomy ID. The resulting file is therefore a concatenated JSON stream rather than one conventional JSON object.

## Stage 2 — Convert concatenated JSON to FASTA

The runner executes `json_to_fasta.py` on `retrieved_proteins.json`.

Output:

```text
converted.fasta
```

The converter uses repeated JSON decoding so multiple concatenated JSON objects can be read. FASTA headers retain accession/title/organism information when available.

## Stage 3 — Prepare alignment input

The supplied `run_alignments.py` expects:

```text
~/Downloads/converted.fasta
```

The integrated runner handles this automatically by temporarily copying the run-specific FASTA to that location.

If a file already exists there, it is backed up before the pipeline copy is installed. After a successful run, the original file is restored; if there was no original file, the temporary copy is removed.

## Stage 4 — MAFFT multiple-sequence alignment

`run_alignments.py` groups sequences using FASTA-header annotations for:

```text
nsp5_3CLpro
nsp12_RdRp
nsp13_helicase
nsp14_ExoN
nsp15_NendoU
nsp16_2O_MTase
```

Groups with fewer than two sequences are skipped.

MAFFT is called with:

```bash
mafft --auto --thread -1 input.fasta
```

Outputs are written under:

```text
alignments/
```

including files such as:

```text
nsp5_3CLpro_aligned.fasta
nsp12_RdRp_aligned.fasta
nsp16_2O_MTase_aligned.fasta
```

At least one `*_aligned.fasta` must be produced for the pipeline to continue.

## Stage 5 — Pairwise identity

`pairwise_identity.py` processes every:

```text
alignments/*_aligned.fasta
```

For every unique sequence pair, it calculates:

```text
percent identity =
    identical residues
    ---------------------- × 100
    compared non-gap sites
```

Positions containing a gap in either sequence are excluded.

Output:

```text
identity_results/all_pairwise_identity_ascending.csv
```

The rows are sorted from lowest to highest pairwise identity.

## Stage 6 — Alignment conservation analysis

`analyze_conservation.py` analyzes every aligned FASTA and calculates per-column consensus, conservation fraction, gap fraction, and complete identity.

Its configured highly conserved criterion is:

```text
conservation >= 0.90
AND
gap_fraction == 0
```

Outputs include:

```text
conservation_results/
├── conservation_summary.csv
├── <protein>_conservation.csv
└── <protein>_conserved_blocks.txt
```

## Stage 7 — Conserved peptide-region extraction

The integrated runner executes `conserved_regions.py` once for every:

```text
alignments/*_aligned.fasta
```

For example:

```bash
python conserved_regions.py \
    alignments/nsp16_2O_MTase_aligned.fasta \
    --min-identity 0.95 \
    --min-occupancy 0.50 \
    --min-length 10 \
    --max-length 15 \
    --weighting none \
    --matrix-out conserved_regions/nsp16_2O_MTase_aligned_matrix.tsv \
    --fasta-out conserved_regions/nsp16_2O_MTase_aligned_conserved_regions.fasta
```

A column is considered conserved when its most common standard amino acid reaches the configured identity threshold and its non-gap occupancy reaches the configured occupancy threshold.

### Matrix output

For each alignment:

```text
<protein>_aligned_matrix.tsv
```

contains sequence-specific information for each conserved stretch, including aligned sequence slices, residue spans, and identity relative to the consensus.

Example:

```text
conserved_regions/nsp16_2O_MTase_aligned_matrix.tsv
```

### Conserved-region FASTA output

For each alignment:

```text
<protein>_aligned_conserved_regions.fasta
```

contains the actual ungapped subsequence carried by each input sequence for each conserved stretch.

FASTA headers include metadata such as:

```text
region=<alignment-span>
residues=<sequence-residue-span>
identity=<identity-to-consensus>
```

Example:

```text
conserved_regions/nsp16_2O_MTase_aligned_conserved_regions.fasta
```

These FASTA peptide sequences become the primary query sequences for the next stage.

## Stage 8 — Search conserved peptides in `cov_pos_ab.tsv`

The final stage searches peptide candidates generated by `conserved_regions.py` against the `peptide` column in:

```text
cov_pos_ab.tsv
```

The database path is controlled by:

```bash
--cov-pos-ab
```

If omitted, the runner expects:

```text
cov_pos_ab.tsv
```

in the project directory.

### Matching rule

The pipeline performs a:

```text
case-insensitive exact peptide sequence match
```

The FASTA sequence is converted to uppercase, and values from the database `peptide` column are stripped and converted to uppercase before lookup.

This is an **exact sequence lookup**. It is not:

- substring matching;
- approximate/fuzzy matching;
- sequence alignment;
- similarity searching; or
- BLAST searching.

For example, a candidate:

```text
ABCDEFGHIJK
```

matches a database peptide only when the complete normalized database value is also:

```text
ABCDEFGHIJK
```

### Input files used by the search

The runner discovers all:

```text
conserved_regions/*_conserved_regions.fasta
```

files.

For each peptide FASTA, it also looks for the corresponding:

```text
conserved_regions/*_matrix.tsv
```

file. The matrix is read to retain the association with the conserved-region analysis, while region/residue/identity metadata used in the match output are parsed from the conserved-region FASTA headers.

### Peptide database requirement

`cov_pos_ab.tsv` must be a tab-separated file containing a column named exactly:

```text
peptide
```

If the file is missing or that column is absent, the pipeline stops with an error.

### Match output

For each conserved-region FASTA, the runner writes:

```text
cov_pos_ab_matches/<protein>_aligned_cov_pos_ab_matches.tsv
```

For example:

```text
cov_pos_ab_matches/nsp16_2O_MTase_aligned_cov_pos_ab_matches.tsv
```

The output fields currently include:

```text
source_fasta
sequence_id
peptide
peptide_length
region
residues
identity
matched
structure_id
source_organism
protein
protein_all
receptor
receptor_confidence
is_bcell
is_tcell
protein_role
protein_receptor
```

For a matching peptide:

```text
matched = 1
```

Database metadata are copied into the corresponding result row.

For a candidate peptide with no exact database match:

```text
matched = 0
```

The candidate is still retained in the result table, while database-derived metadata fields remain empty.

If one peptide occurs in multiple rows of `cov_pos_ab.tsv`, the output can contain multiple rows for that candidate—one for each matching database row.

### Match summary

The runner also writes:

```text
cov_pos_ab_matches/cov_pos_ab_match_summary.tsv
```

with:

```text
source_fasta
candidates
matched_candidates
output
```

This provides a per-conserved-region-FASTA summary and links each source to its detailed match table.

## Output directory structure

A successful run resembles:

```text
pipeline_runs/
└── YYYYMMDD_HHMMSS/
    ├── retrieved_proteins.json
    ├── converted.fasta
    ├── pipeline.log
    ├── pipeline_manifest.json
    │
    ├── alignments/
    │   ├── nsp5_3CLpro.fasta
    │   ├── nsp5_3CLpro_aligned.fasta
    │   ├── nsp16_2O_MTase.fasta
    │   └── nsp16_2O_MTase_aligned.fasta
    │
    ├── identity_results/
    │   └── all_pairwise_identity_ascending.csv
    │
    ├── conservation_results/
    │   ├── conservation_summary.csv
    │   ├── *_conservation.csv
    │   └── *_conserved_blocks.txt
    │
    ├── conserved_regions/
    │   ├── nsp5_3CLpro_aligned_matrix.tsv
    │   ├── nsp5_3CLpro_aligned_conserved_regions.fasta
    │   ├── nsp16_2O_MTase_aligned_matrix.tsv
    │   └── nsp16_2O_MTase_aligned_conserved_regions.fasta
    │
    └── cov_pos_ab_matches/
        ├── nsp5_3CLpro_aligned_cov_pos_ab_matches.tsv
        ├── nsp16_2O_MTase_aligned_cov_pos_ab_matches.tsv
        └── cov_pos_ab_match_summary.tsv
```

## Pipeline manifest

A successful run creates:

```text
pipeline_manifest.json
```

The manifest records run settings and output locations, including the NCBI JSON, converted FASTA, alignment files, pairwise identity results, conservation results, conserved-region matrix/FASTA outputs, and `cov_pos_ab.tsv` peptide-search results.

This makes it easier to identify all files produced by one run.

## Pipeline log

Subprocess output is streamed to the terminal and recorded in:

```text
pipeline.log
```

The log records commands and exit codes. A nonzero subprocess exit status stops the pipeline.

## Recommended project layout

```text
viral-analysis/
├── integrated_viral_pipeline.py
├── viral_protein_pipeline.py
├── json_to_fasta.py
├── run_alignments.py
├── pairwise_identity.py
├── analyze_conservation.py
├── conserved_regions.py
├── cov_pos_ab.tsv
├── taxonomy_ids.csv
└── requirements.txt
```

Run:

```bash
cd viral-analysis

python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --cov-pos-ab cov_pos_ab.tsv
```

## Conda setup example

```bash
conda create -n viral-pipeline python=3.11 pandas numpy biopython
conda activate viral-pipeline
conda install -c bioconda mafft
```

Then:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --cov-pos-ab cov_pos_ab.tsv
```

## Pre-run checks

```bash
python --version
python -c "import pandas, numpy, Bio; print('Python dependencies OK')"
mafft --version
```

Confirm these files are present:

```text
integrated_viral_pipeline.py
viral_protein_pipeline.py
json_to_fasta.py
run_alignments.py
pairwise_identity.py
analyze_conservation.py
conserved_regions.py
cov_pos_ab.tsv
taxonomy_ids.csv
```

Confirm the Taxonomy ID file contains:

```text
Taxonomy ID
```

and `cov_pos_ab.tsv` contains:

```text
peptide
```

## Troubleshooting

### `Required script not found`

The integrated runner looks for the component scripts in its own directory using their exact filenames.

### `ModuleNotFoundError: No module named 'pandas'`

```bash
python -m pip install pandas
```

### `ModuleNotFoundError: No module named 'numpy'`

```bash
python -m pip install numpy
```

### `ModuleNotFoundError: No module named 'Bio'`

```bash
python -m pip install biopython
```

### `mafft` not found

Verify:

```bash
mafft --version
```

and:

```bash
python -c "import shutil; print(shutil.which('mafft'))"
```

### Taxonomy CSV error

The input CSV must contain a column named exactly:

```text
Taxonomy ID
```

### `cov_pos_ab.tsv not found`

Provide the path explicitly:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --cov-pos-ab /full/path/to/cov_pos_ab.tsv
```

### `does not contain a 'peptide' column`

The peptide database must be tab-separated and contain a column named exactly:

```text
peptide
```

### No alignment files produced

`run_alignments.py` skips protein groups containing fewer than two matching sequences. If no target group has at least two sequences, no `*_aligned.fasta` files are generated and the integrated runner stops.

### No conserved-region FASTA files

The peptide lookup stage expects files matching:

```text
conserved_regions/*_conserved_regions.fasta
```

These are generated by the conserved-region stage.

A FASTA file can exist but contain no peptide records if the alignment does not contain stretches satisfying the configured conservation and minimum-length criteria. In that situation, its detailed peptide match output will contain no candidate records.

### No matches in `cov_pos_ab.tsv`

An unmatched peptide is not necessarily an error. The lookup is exact. Differences in length or even one amino acid prevent a match.

Review the detailed output and check:

```text
matched
```

where `0` indicates no exact peptide match.

### Existing work directory

Use a new directory or explicitly permit reuse:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --work-dir results/run1 \
    --keep-existing-run-dir
```

## Interpretation notes

Pairwise identity excludes positions containing a gap in either sequence.

`analyze_conservation.py` uses its own 90%/zero-gap conservation definition.

`conserved_regions.py` uses the separately configurable conserved-column identity and occupancy thresholds and extracts peptide-length windows/stretches from those regions.

The final `cov_pos_ab.tsv` step performs only exact normalized sequence matching against the `peptide` column. A database match indicates that the same peptide sequence is present in the supplied database; it does not by itself establish functional or experimental equivalence beyond the metadata present in that database.

## Quick start

```bash
# Install Python dependencies
python -m pip install pandas numpy biopython

# Install/verify MAFFT
mafft --version

# Confirm taxonomy_ids.csv has a "Taxonomy ID" column
# Confirm cov_pos_ab.tsv has a "peptide" column

python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --cov-pos-ab cov_pos_ab.tsv \
    --limit 200

# Inspect the new timestamped directory under pipeline_runs/
# In particular:
#   conserved_regions/
#   cov_pos_ab_matches/
#   cov_pos_ab_matches/cov_pos_ab_match_summary.tsv
```
