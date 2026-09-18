# Integrated Viral Protein Analysis Pipeline

`integrated_viral_pipeline.py` runs an end-to-end viral protein sequence
workflow that connects sequence retrieval from NCBI with FASTA
conversion, MAFFT multiple-sequence alignment, pairwise identity
analysis, conservation analysis, and conserved peptide-window extraction.

## Workflow

The integrated pipeline executes the following stages in order:

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
run_alignments.py
      |
      v
alignments/*_aligned.fasta
      |
      +---------------------------+------------------------------+
      |                           |                              |
      v                           v                              v
pairwise_identity.py      analyze_conservation.py       conserved_regions.py
      |                           |                    (one run per alignment)
      v                           v                              |
identity_results/         conservation_results/                v
                                                       conserved_regions/
                                                       matrix TSV + FASTA
```

The `pairwise_identity.py`, `analyze_conservation.py`, and `conserved_regions.py` stages all consume the aligned FASTA files produced by `run_alignments.py`. The integrated runner executes `conserved_regions.py` once for every file matching `alignments/*_aligned.fasta`.

The integrated runner stops if a required stage fails or if an expected output file is not produced.

# Pipeline components

The following files must be in the **same directory** as
`integrated_viral_pipeline.py`:

``` text
integrated_viral_pipeline.py
viral_protein_pipeline.py
json_to_fasta.py
run_alignments.py
pairwise_identity.py
analyze_conservation.py
conserved_regions.py
```

The filenames are significant because the integrated runner currently
looks for these exact names.

## Requirements

### Python

Python **3.9 or newer** is recommended.

The integrated runner and supporting scripts use modern Python features
such as built-in generic type annotations:

``` python
list[str]
```

Check your Python version with:

``` bash
python --version
```

or:

``` bash
python3 --version
```

### Python packages

The sequence-retrieval script requires `pandas` to read the Taxonomy ID
CSV.

Install it with:

``` bash
python -m pip install pandas
```

The other supplied downstream scripts use Python standard-library
modules only.

A minimal `requirements.txt` for this integrated workflow can therefore
contain:

``` text
pandas>=2.0
```

Install from it with:

``` bash
python -m pip install -r requirements.txt
```

### MAFFT

`run_alignments.py` calls the external `mafft` executable, so MAFFT must
be installed separately and available on your system `PATH`.

Verify installation with:

``` bash
mafft --version
```

If the command is not found, install MAFFT.

#### macOS with Homebrew

``` bash
brew install mafft
```

#### Ubuntu/Debian

``` bash
sudo apt update
sudo apt install mafft
```

#### Conda/Bioconda

``` bash
conda install -c bioconda mafft
```

or:

``` bash
mamba install -c bioconda mafft
```

You can also verify that Python can locate MAFFT:

``` bash
python -c "import shutil; print(shutil.which('mafft'))"
```

A valid installation should print a path to the executable rather than
`None`.

### Conserved-region analysis dependencies

`conserved_regions.py` requires:

- **NumPy** (`numpy`)
- **Biopython** (`biopython`)

Install them with:

```bash
python -m pip install numpy biopython
```

A minimal Python dependency set for the integrated sequence workflow is therefore:

```text
pandas>=2.0
numpy
biopython
```

`MAFFT` remains an external command-line dependency and is not installed by `pip`.

### Internet access

Internet access is required during the fetch stage because
`viral_protein_pipeline.py` queries NCBI E-utilities.

The retrieval code uses NCBI Protein ESearch, ESummary, and EFetch
endpoints.

## Input Taxonomy ID CSV

The required pipeline input is a CSV file containing a column named
exactly:

``` text
Taxonomy ID
```

Example `taxonomy_ids.csv`:

``` csv
Taxonomy ID
2697049
694009
11137
```

Multiple Taxonomy IDs may be supplied.

The fetch implementation reads the CSV with pandas and executes the NCBI
retrieval step for every value in the `Taxonomy ID` column.

## Basic usage

Run the pipeline with:

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv
```

On systems where the interpreter is named `python3`:

``` bash
python3 integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv
```

## Command-line arguments

  ---------------------------------------------------------------------------------
  Argument                    Required          Default           Description
  --------------------------- ----------------- ----------------- -----------------
  `--taxid-csv`               Yes               ---               CSV containing a
                                                                  `Taxonomy ID`
                                                                  column.

  `--limit`                   No                `200`             Maximum NCBI
                                                                  records requested
                                                                  per Taxonomy ID.

  `--host-filter`             No                Off               Enables the host
                                                                  filter
                                                                  implemented by
                                                                  the retrieval
                                                                  script.

  `--work-dir`                No                Timestamped       Directory used
                                                directory         for all outputs
                                                                  from the run.

  `--python`                  No                Current           Python executable
                                                interpreter       used to execute
                                                                  the component
                                                                  scripts.

  `--fasta-width`             No                `60`              Number of
                                                                  residues per
                                                                  FASTA sequence
                                                                  line.

  `--keep-existing-run-dir`   No                Off               Allows reuse of
                                                                  an already
                                                                  existing work
                                                                  directory.
  ---------------------------------------------------------------------------------

Display the built-in help:

``` bash
python integrated_viral_pipeline.py --help
```

## Example: increase retrieval limit

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --limit 500
```

This requests up to 500 protein records for each Taxonomy ID.

## Example: enable host filtering

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --limit 500 \
    --host-filter
```

## Example: choose an output directory

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --work-dir results/coronavirus_run
```

If the directory already exists, the runner stops unless:

``` bash
--keep-existing-run-dir
```

is supplied.

For example:

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --work-dir results/coronavirus_run \
    --keep-existing-run-dir
```

## Stage 1 --- Retrieve proteins from NCBI

The integrated runner executes the equivalent of:

``` bash
python "viral_protein_pipeline(1).py" fetch \
    --taxid taxonomy_ids.csv \
    --limit 200 \
    --output retrieved_proteins.json
```

For every Taxonomy ID, the retrieval script searches the NCBI Protein
database for RefSeq protein records and retrieves sequence and metadata
fields including:

``` text
uid
accession
title
length
organism
sequence
```

The output is:

``` text
retrieved_proteins.json
```

### Concatenated JSON behavior

The retrieval function opens the output file in append mode. Therefore,
when several Taxonomy IDs are present, the file contains multiple JSON
objects one after another.

Conceptually:

``` text
{ JSON object for TaxID 1 }
{ JSON object for TaxID 2 }
{ JSON object for TaxID 3 }
```

This is why the next stage uses the supplied `json_to_fasta.py`, which
is designed to parse concatenated JSON objects.

The integrated runner also counts the JSON objects after retrieval and
reports the count in the pipeline log.

## Stage 2 --- Convert JSON to FASTA

The runner next executes:

``` bash
python json_to_fasta.py \
    retrieved_proteins.json \
    --output converted.fasta \
    --width 60
```

`json_to_fasta.py` repeatedly uses `JSONDecoder.raw_decode()` so that
concatenated JSON objects can be processed.

For each valid record, the FASTA header contains the accession, title,
and organism when available.

Conceptually:

``` text
>ACCESSION protein title | organism
AMINOACIDSEQUENCE...
```

The default sequence line width is 60 residues.

Change it with:

``` bash
--fasta-width 80
```

at the integrated-pipeline level.

## Stage 3 --- Prepare the alignment input

The supplied `run_alignments.py` currently expects its input FASTA at:

``` text
~/Downloads/converted.fasta
```

The integrated runner handles this automatically.

Before alignment it:

1.  creates `~/Downloads` if necessary;
2.  backs up an existing `~/Downloads/converted.fasta`;
3.  copies the run-specific `converted.fasta` to that location;
4.  executes `run_alignments.py`; and
5.  after a successful pipeline, restores the previous file or removes
    the temporary copy.

The run-specific FASTA remains safely stored inside the pipeline output
directory.

## Stage 4 --- Multiple-sequence alignment with MAFFT

`run_alignments.py` reads the combined FASTA and groups proteins
according to patterns in their FASTA headers.

The configured target groups are:

``` text
nsp5_3CLpro
nsp12_RdRp
nsp13_helicase
nsp14_ExoN
nsp15_NendoU
nsp16_2O_MTase
```

A group containing fewer than two sequences is skipped.

For eligible groups, MAFFT is executed approximately as:

``` bash
mafft --auto --thread -1 input.fasta
```

`--auto` allows MAFFT to select an alignment strategy, while
`--thread -1` lets MAFFT automatically use available CPU threads.

Outputs are stored under:

``` text
alignments/
```

For example:

``` text
alignments/
├── nsp5_3CLpro.fasta
├── nsp5_3CLpro_aligned.fasta
├── nsp12_RdRp.fasta
├── nsp12_RdRp_aligned.fasta
└── ...
```

The integrated runner requires at least one:

``` text
*_aligned.fasta
```

file before proceeding.

## Stage 5 --- Pairwise sequence identity

`pairwise_identity.py` reads:

``` text
alignments/*_aligned.fasta
```

and calculates amino-acid identity for every unique sequence pair within
each protein alignment.

Alignment columns containing a gap in either sequence are excluded from
the identity denominator.

The calculation is:

``` text
percent identity =
    identical residues
    -------------------  x 100
    compared non-gap positions
```

The results are sorted from lowest to highest identity and written to:

``` text
identity_results/all_pairwise_identity_ascending.csv
```

Columns include:

``` text
protein
sequence_1
sequence_2
identical_residues
positions_compared
percent_identity
```

## Stage 6 --- Conservation analysis

`analyze_conservation.py` also reads:

``` text
alignments/*_aligned.fasta
```

For each alignment column it calculates:

-   consensus residue;
-   conservation fraction;
-   gap fraction; and
-   whether the position is completely identical across all sequences.

The configured highly conserved threshold is:

``` text
0.90
```

A position is considered highly conserved when:

``` text
conservation >= 0.90
AND
gap_fraction == 0
```

Consecutive conserved positions are reported as conserved blocks when
the block contains at least three positions.

### Conservation outputs

The script creates:

``` text
conservation_results/
```

For each protein, it writes a position-level CSV:

``` text
<protein>_conservation.csv
```

with columns:

``` text
alignment_position
consensus
conservation
gap_fraction
identical
```

It also writes:

``` text
<protein>_conserved_blocks.txt
```

containing the coordinates, lengths, and consensus sequences of
conserved blocks.

A combined summary is written to:

``` text
conservation_results/conservation_summary.csv
```

with columns:

``` text
protein
sequences
alignment_length
identical_sites
highly_conserved_sites
percent_identical_sites
percent_highly_conserved
conserved_blocks
```

## Output directory structure

Unless `--work-dir` is specified, the integrated runner creates:

``` text
pipeline_runs/YYYYMMDD_HHMMSS/
```

A successful run will resemble:

``` text
pipeline_runs/
└── 20260917_145100/
    ├── retrieved_proteins.json
    ├── converted.fasta
    ├── pipeline.log
    ├── pipeline_manifest.json
    │
    ├── alignments/
    │   ├── nsp5_3CLpro.fasta
    │   ├── nsp5_3CLpro_aligned.fasta
    │   ├── nsp12_RdRp.fasta
    │   ├── nsp12_RdRp_aligned.fasta
    │   └── ...
    │
    ├── identity_results/
    │   └── all_pairwise_identity_ascending.csv
    │
    └── conservation_results/
        ├── conservation_summary.csv
        ├── nsp5_3CLpro_conservation.csv
        ├── nsp5_3CLpro_conserved_blocks.txt
        └── ...
```

## Pipeline manifest

After a successful run, the runner creates:

``` text
pipeline_manifest.json
```

The manifest records the principal run settings and output locations,
including:

``` text
taxid_csv
limit
host_filter
fasta_width
work_dir
JSON output
FASTA output
alignment files
pairwise identity output
conservation summary
conservation output directory
```

This file is useful for tracking the products of a particular run.

## Pipeline log

All subprocess output is streamed to the terminal and also written to:

``` text
pipeline.log
```

The log records each executed command and its exit code.

If a stage returns a nonzero exit status, the integrated runner stops
rather than continuing with incomplete data.

## Recommended project layout

Place all scripts together:

``` text
viral-analysis/
├── integrated_viral_pipeline.py
├── viral_protein_pipeline(1).py
├── json_to_fasta.py
├── run_alignments.py
├── pairwise_identity.py
├── analyze_conservation.py
├── requirements.txt
└── taxonomy_ids.csv
```

Then:

``` bash
cd viral-analysis
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv
```

## Suggested `requirements.txt`

For the supplied scripts, a minimal Python requirements file is:

``` text
pandas>=2.0
```

MAFFT should **not** be placed in `requirements.txt`, because it is an
external command-line application rather than a Python package.

Install the two dependency layers separately:

``` bash
python -m pip install -r requirements.txt
```

and, for example on macOS:

``` bash
brew install mafft
```

## Conda environment example

A convenient bioinformatics environment can be created with Conda:

``` bash
conda create -n viral-pipeline python=3.11 pandas
conda activate viral-pipeline
conda install -c bioconda mafft
```

Then run:

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv
```

## Pre-run checks

Before starting a large run, verify:

``` bash
python --version
python -c "import pandas; print(pandas.__version__)"
mafft --version
```

Then verify that all required scripts are present in the same directory:

``` text
integrated_viral_pipeline.py
viral_protein_pipeline(1).py
json_to_fasta.py
run_alignments.py
pairwise_identity.py
analyze_conservation.py
conserved_regions.py
```

Also confirm that the input CSV contains:

``` text
Taxonomy ID
```

as its column header.

## Troubleshooting

### `Required script not found`

The integrated runner searches its own directory for the five component
scripts using their exact filenames.

Make sure the files are named:

``` text
viral_protein_pipeline(1).py
json_to_fasta.py
run_alignments.py
pairwise_identity.py
analyze_conservation.py
conserved_regions.py
```

### `ModuleNotFoundError: No module named 'pandas'`

Install pandas:

``` bash
python -m pip install pandas
```

### `FileNotFoundError: mafft`

or:

``` text
No such file or directory: 'mafft'
```

MAFFT is either not installed or is not on the system `PATH`.

Check:

``` bash
mafft --version
```

and:

``` bash
python -c "import shutil; print(shutil.which('mafft'))"
```

### Taxonomy CSV error

The fetch script expects a CSV column named exactly:

``` text
Taxonomy ID
```

A valid example is:

``` csv
Taxonomy ID
2697049
```

### No alignment files were produced

The alignment script only aligns a target protein group when at least
two matching sequences are available.

The target groups are nsp5, nsp12, nsp13, nsp14, nsp15, and
nsp16-related proteins. If the fetched data contain fewer than two
matching sequences for every target, MAFFT alignment files will not be
generated and the integrated runner will stop.

### Existing work directory

By default, the runner will not reuse an existing `--work-dir`.

Either choose a new directory:

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --work-dir results/run2
```

or explicitly allow reuse:

``` bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --work-dir results/run1 \
    --keep-existing-run-dir
```

### NCBI request failure

The fetch stage requires network access to NCBI E-utilities. Check
internet access, proxy/firewall configuration, and NCBI service
availability.

The retrieval script retries failed HTTP requests before reporting an
error.

## Stage 7 — Conserved peptide-region extraction

After the aligned FASTA files have been created, the runner executes `conserved_regions.py` **once for every** alignment matching:

```text
alignments/*_aligned.fasta
```

For example, a file such as:

```text
alignments/nsp5_3CLpro_aligned.fasta
```

is processed with a command equivalent to:

```bash
python conserved_regions.py alignments/nsp5_3CLpro_aligned.fasta \
    --min-identity 0.95 \
    --min-occupancy 0.5 \
    --min-length 10 \
    --max-length 15 \
    --weighting none \
    --matrix-out conserved_regions/nsp5_3CLpro_aligned_matrix.tsv \
    --fasta-out conserved_regions/nsp5_3CLpro_aligned_conserved_regions.fasta
```

The integrated runner derives the output names from the alignment filename and repeats the process for every other `*_aligned.fasta` file.

### Conserved-region options

| Integrated option | Default | Meaning |
| --- | ---: | --- |
| `--conserved-min-identity` | `0.95` | Minimum fraction of a conserved column sharing one amino acid. Values such as `95` are also accepted. |
| `--conserved-min-occupancy` | `0.50` | Minimum fraction of sequences that must be non-gap at the column. |
| `--conserved-min-length` | `10` | Shortest peptide window reported. |
| `--conserved-max-length` | `15` | Longest peptide window reported. |
| `--conserved-weighting` | `none` | Sequence weighting method: `none` or `henikoff`. |
| `--conserved-gap-votes` | Off | Includes gaps in the identity denominator when enabled. |

For example:

```bash
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --conserved-min-identity 90 \
    --conserved-min-occupancy 0.8 \
    --conserved-min-length 10 \
    --conserved-max-length 20 \
    --conserved-weighting henikoff
```

The underlying `conserved_regions.py` script defines a column as conserved when the most common amino acid meets the identity threshold and column occupancy meets the occupancy threshold. Consecutive conserved columns are then used to generate every peptide window between the configured minimum and maximum lengths.

### `--matrix-out` output

For each alignment, `--matrix-out` creates a tab-separated matrix with one row per input sequence and three fields for each conserved stretch: the aligned sequence slice, the sequence's own residue-number span, and identity relative to the consensus. The consensus motif is also included in the output.

Example:

```text
conserved_regions/nsp5_3CLpro_aligned_matrix.tsv
```

### `--fasta-out` output

For each alignment, `--fasta-out` creates a FASTA file containing the **actual ungapped sequence subsequences** carried by each sequence for each conserved stretch. These records are grouped by conserved stretch and are the appropriate output to inspect when selecting sequence-specific peptide candidates.

Example:

```text
conserved_regions/nsp5_3CLpro_aligned_conserved_regions.fasta
```

The consensus window reported by the window table does not necessarily occur exactly in any one sequence; the matrix and FASTA outputs preserve what individual sequences actually contain.

### Conserved-region output directory

The integrated runner creates:

```text
conserved_regions/
├── nsp5_3CLpro_aligned_matrix.tsv
├── nsp5_3CLpro_aligned_conserved_regions.fasta
├── nsp12_RdRp_aligned_matrix.tsv
├── nsp12_RdRp_aligned_conserved_regions.fasta
├── nsp13_helicase_aligned_matrix.tsv
├── nsp13_helicase_aligned_conserved_regions.fasta
└── ...
```

## Scientific interpretation notes

The downstream analyses operate on protein groups identified from
FASTA-header annotations.

Pairwise identity excludes alignment positions containing a gap in
either member of a pair.

The conservation analysis defines highly conserved positions using a
fixed 90% threshold with no gaps and reports conserved blocks of at
least three consecutive positions.

These definitions are implementation choices of the supplied scripts and
should be considered when interpreting downstream results.

## Quick start

``` bash
# 1. Install Python dependencies
python -m pip install pandas numpy biopython

# 2. Verify MAFFT
mafft --version

# 3. Prepare taxonomy_ids.csv with a "Taxonomy ID" column

# 4. Run
python integrated_viral_pipeline.py \
    --taxid-csv taxonomy_ids.csv \
    --limit 200

# 5. Inspect the timestamped directory under pipeline_runs/
```
