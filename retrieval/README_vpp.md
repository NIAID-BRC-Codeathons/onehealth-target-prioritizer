# Viral Protein Processing Pipeline

`viral_protein_pipeline.py` is a command-line Python pipeline for
retrieving viral protein sequences from NCBI, performing sequence
quality control, classifying structural proteins, and identifying
conserved regions across protein sequences.

The pipeline has four commands:

1.  `fetch` --- retrieve RefSeq protein records from NCBI.
2.  `qc` --- filter sequences by length, ambiguity, and duplication.
3.  `annotate` --- classify proteins into S, M, N, E, and unassigned
    groups.
4.  `msa-conservation` --- perform a built-in progressive alignment and
    identify conserved sequence blocks.

## Requirements

-   Python 3
-   `pandas`
-   Internet access for the `fetch` command

Install the Python dependency with:

``` bash
python -m pip install pandas
```

All other imported modules are from the Python standard library.

## Important input convention

Despite the `--taxid` option name, the current `fetch` implementation
expects `--taxid` to point to a **CSV file**, not to a single numeric
Taxonomy ID.

The CSV must contain a column named:

``` text
Taxonomy ID
```

For example:

``` csv
Taxonomy ID
2697049
694009
11137
```

The script reads that column with pandas and runs the NCBI retrieval
step once for every Taxonomy ID.

## Basic help

Display the available commands:

``` bash
python viral_protein_pipeline.py --help
```

Display help for a particular command:

``` bash
python viral_protein_pipeline.py fetch --help
python viral_protein_pipeline.py qc --help
python viral_protein_pipeline.py annotate --help
python viral_protein_pipeline.py msa-conservation --help
```

## 1. Fetch protein sequences from NCBI

### Usage

``` bash
python viral_protein_pipeline.py fetch \
    --taxid taxonomy_ids.csv \
    --limit 200 \
    --output retrieved_proteins.json
```

With the optional host filter:

``` bash
python viral_protein_pipeline.py fetch \
    --taxid taxonomy_ids.csv \
    --limit 200 \
    --host-filter \
    --output retrieved_proteins.json
```

### Arguments

  ----------------------------------------------------------------------------------
  Argument          Required          Default           Description
  ----------------- ----------------- ----------------- ----------------------------
  `--taxid`         Yes               ---               CSV file containing a
                                                        `Taxonomy ID` column.

  `--limit`         No                `200`             Maximum number of records
                                                        requested for each Taxonomy
                                                        ID.

  `--host-filter`   No                Off               Adds
                                                        mammalian/pangolin-related
                                                        host filtering to the NCBI
                                                        query.

  `--output`        Yes               ---               JSON output path.
  ----------------------------------------------------------------------------------

The retrieval step queries the NCBI Protein database through Entrez
E-utilities and restricts the search to RefSeq protein records.

For each retrieved protein, the script stores fields including:

``` json
{
  "uid": "...",
  "accession": "...",
  "title": "...",
  "length": 123,
  "organism": "...",
  "sequence": "..."
}
```

### Important behavior with multiple Taxonomy IDs

`fetch_viral_proteins()` opens the output file in **append mode**.
Therefore, when the CSV contains multiple Taxonomy IDs, multiple
complete JSON objects are written consecutively to the same output file.

Conceptually, the result is:

``` text
{ JSON object for TaxID 1 }
{ JSON object for TaxID 2 }
{ JSON object for TaxID 3 }
```

This is a concatenated JSON stream rather than a single standard JSON
document.

Consequently, the pipeline's current `qc` command, which uses
`json.load()`, cannot directly read a multi-TaxID fetch output and may
fail with:

``` text
json.decoder.JSONDecodeError: Extra data
```

For the current implementation, use a single Taxonomy ID per fetch
output when the output will be passed directly to `qc`, or
preprocess/merge the concatenated objects into a valid JSON document
before running QC.

Also note that rerunning `fetch` with an existing output filename
appends additional JSON content rather than replacing the existing file.
Remove or rename an old output file before a fresh retrieval if you do
not want that behavior.

## 2. Quality control

The QC step reads a retrieved JSON object and filters its protein
records.

### Usage

``` bash
python viral_protein_pipeline.py qc \
    --input retrieved_proteins.json \
    --output qc_proteins.json
```

Custom sequence-length limits:

``` bash
python viral_protein_pipeline.py qc \
    --input retrieved_proteins.json \
    --output qc_proteins.json \
    --min-length 50 \
    --max-length 3000
```

### Arguments

  Argument         Required   Default   Description
  ---------------- ---------- --------- -----------------------------------------------
  `--input`        Yes        ---       Input JSON file containing a `records` array.
  `--output`       Yes        ---       Output JSON file for sequences passing QC.
  `--min-length`   No         `50`      Minimum accepted sequence length.
  `--max-length`   No         `3000`    Maximum accepted sequence length.

### QC criteria

A sequence is retained when it satisfies the configured length range and
the ambiguity threshold, and it is not an exact sequence duplicate.

The code treats residues outside:

``` text
ACDEFGHIKLMNPQRSTVWY
```

as ambiguous. This includes characters such as `X`, `Z`, `B`, and `*`.

The underlying QC function has a default maximum ambiguous-residue
percentage of `1.0%`. The current command-line interface does not expose
this threshold as an argument.

Exact duplicate sequences are removed using the sequence itself as the
deduplication key.

### QC output

The output contains QC statistics plus the retained records:

``` json
{
  "qc_stats": {
    "total": 100,
    "failed_length": 10,
    "failed_ambiguity": 2,
    "failed_duplicate": 5,
    "passed": 83
  },
  "records": []
}
```

## 3. Annotate and classify proteins

The annotation step classifies records according to regular-expression
matches against each record's `title`.

### Usage

``` bash
python viral_protein_pipeline.py annotate \
    --input qc_proteins.json \
    --output-dir annotated
```

### Protein classes

The script attempts to assign proteins to:

  Class          Intended category
  -------------- --------------------------------------------
  `S`            Spike / surface glycoprotein
  `M`            Membrane / matrix protein
  `N`            Nucleocapsid / nucleoprotein
  `E`            Envelope / small membrane protein
  `unassigned`   Titles that do not match the S/M/N/E rules

Classification is title-based rather than sequence-based.

### Output files

For each class, the script writes both JSON and FASTA files:

``` text
annotated/
├── proteins_S.json
├── proteins_S.fasta
├── proteins_M.json
├── proteins_M.fasta
├── proteins_N.json
├── proteins_N.fasta
├── proteins_E.json
├── proteins_E.fasta
├── proteins_unassigned.json
└── proteins_unassigned.fasta
```

A FASTA entry has the form:

``` text
>ACCESSION protein title
SEQUENCE
```

## 4. Multiple sequence alignment and conserved-region extraction

The final stage processes the S, M, N, and E JSON files produced by
`annotate`.

### Usage

``` bash
python viral_protein_pipeline.py msa-conservation \
    --input-dir annotated \
    --output-dir conservation
```

Custom conservation parameters:

``` bash
python viral_protein_pipeline.py msa-conservation \
    --input-dir annotated \
    --output-dir conservation \
    --min-region-len 15 \
    --max-entropy 0.3
```

### Arguments

  -----------------------------------------------------------------------------
  Argument             Required          Default           Description
  -------------------- ----------------- ----------------- --------------------
  `--input-dir`        Yes               ---               Directory containing
                                                           `proteins_S.json`,
                                                           `proteins_M.json`,
                                                           `proteins_N.json`,
                                                           and
                                                           `proteins_E.json`.

  `--output-dir`       Yes               ---               Directory for
                                                           conservation
                                                           results.

  `--min-region-len`   No                `15`              Minimum number of
                                                           consecutive
                                                           conserved alignment
                                                           columns required for
                                                           a region.

  `--max-entropy`      No                `0.3`             Maximum Shannon
                                                           entropy allowed for
                                                           a conserved column.
  -----------------------------------------------------------------------------

Classes containing fewer than two records are skipped.

## Alignment method

The script implements alignment internally and does not call an external
program such as MAFFT, MUSCLE, or Clustal Omega.

It uses:

-   Needleman-Wunsch pairwise global alignment;
-   match score: `+2`;
-   mismatch score: `-1`;
-   gap score: `-2`; and
-   a sequential progressive procedure anchored on the first sequence.

Because this is a simple built-in implementation, computational cost can
become substantial for numerous or long viral proteins.

## Conservation definition

For every alignment column, the pipeline calculates:

-   Shannon entropy;
-   gap ratio; and
-   consensus residue.

A column is marked conserved when:

``` text
entropy <= max_entropy
AND
gap_ratio <= 0.20
```

Consecutive conserved columns are combined into a conserved region when
the block contains at least `--min-region-len` positions.

For each retained region, the report records:

-   alignment start position;
-   alignment end position;
-   region length;
-   consensus motif; and
-   average entropy.

## Conservation output

The final report is written to:

``` text
conservation/conserved_regions.json
```

Its structure is approximately:

``` json
{
  "S": {
    "num_sequences": 10,
    "alignment_length": 1400,
    "conserved_regions_count": 2,
    "conserved_regions": [
      {
        "start": 100,
        "end": 120,
        "length": 21,
        "motif": "EXAMPLECONSENSUSMOTIF",
        "avg_entropy": 0.12
      }
    ]
  }
}
```

Positions reported by the conservation step are **alignment positions**,
not necessarily residue coordinates in an individual ungapped protein.

## Complete example workflow

Create a Taxonomy ID CSV:

``` csv
Taxonomy ID
2697049
```

Then run:

``` bash
# Step 1: Retrieve proteins
python viral_protein_pipeline.py fetch \
    --taxid taxonomy_ids.csv \
    --limit 200 \
    --output retrieved_proteins.json

# Step 2: Quality control
python viral_protein_pipeline.py qc \
    --input retrieved_proteins.json \
    --output qc_proteins.json \
    --min-length 50 \
    --max-length 3000

# Step 3: Structural-protein classification
python viral_protein_pipeline.py annotate \
    --input qc_proteins.json \
    --output-dir annotated

# Step 4: Alignment and conservation analysis
python viral_protein_pipeline.py msa-conservation \
    --input-dir annotated \
    --output-dir conservation \
    --min-region-len 15 \
    --max-entropy 0.3
```

The resulting directory might look like:

``` text
project/
├── viral_protein_pipeline.py
├── taxonomy_ids.csv
├── retrieved_proteins.json
├── qc_proteins.json
├── annotated/
│   ├── proteins_S.json
│   ├── proteins_S.fasta
│   ├── proteins_M.json
│   ├── proteins_M.fasta
│   ├── proteins_N.json
│   ├── proteins_N.fasta
│   ├── proteins_E.json
│   ├── proteins_E.fasta
│   ├── proteins_unassigned.json
│   └── proteins_unassigned.fasta
└── conservation/
    └── conserved_regions.json
```

## NCBI access

The retrieval stage uses the NCBI E-utilities endpoints for:

-   ESearch;
-   ESummary; and
-   EFetch.

HTTP requests use a custom user-agent and retry failed requests up to
three times with exponential backoff.

The script does not currently provide command-line options for an NCBI
API key or email address.

## Known implementation considerations

1.  **`--taxid` is a CSV path.** The CLI help text describes it as an
    NCBI Taxonomy ID, but the implementation passes the value to
    `pandas.read_csv()` and expects a `Taxonomy ID` column.

2.  **Multi-TaxID retrieval creates concatenated JSON.** Each Taxonomy
    ID appends a separate JSON object to the same output file. This
    output is not directly compatible with the current `qc`
    implementation, which expects one JSON object.

3.  **Fetch output is appended.** Existing output files are not
    automatically cleared.

4.  **QC ambiguity threshold is not exposed through the CLI.** The
    function supports `max_ambiguous_pct`, but the command-line parser
    does not currently define an option for it.

5.  **Annotation is title-based.** Classification depends on
    regular-expression matches against NCBI protein titles.

6.  **The MSA implementation is intentionally simple.** It is a custom
    progressive procedure based on pairwise Needleman-Wunsch alignment
    and should not be assumed to reproduce established MSA software.

## License

No license information is specified in the supplied script. Add an
appropriate license section here if the project is distributed or
published.

## Notes
- README generated by GPT-5.6 SoI
- Script authored by Milana Djonovic, co-authored by AntiGravity IDE, modified by Sachin Kumar (to accept file with multiple tax ids)