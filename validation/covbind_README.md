# CovBind Literature Pipeline

`covbind_pipeline_openai.py` is a literature-mining pipeline for finding
papers about **protein binding in alpha- and betacoronaviruses**,
optionally filtering papers by journal metrics, extracting structured
binding information from abstracts, and exporting the results to CSV and
Excel.

## What the pipeline does

The script performs four main tasks:

1.  Searches **Europe PMC** for alpha- and betacoronavirus
    protein-binding literature.
2.  Optionally filters papers using journal metrics from
    `scimagojr.csv`.
3.  Extracts the viral protein, protein region/domain, binding partner,
    and biological importance from each abstract using the configured
    OpenAI-compatible API.
4.  Writes the resulting dataset to `covbind_results.csv` and
    `covbind_results.xlsx`.

If model-based extraction fails, the script falls back to a limited
keyword and regular-expression extractor.

## Project files

A typical project directory is:

``` text
project/
├── covbind_pipeline_openai.py
├── requirements.txt
└── scimagojr.csv          # optional
```

After a successful run:

``` text
project/
├── covbind_pipeline_openai.py
├── requirements.txt
├── scimagojr.csv
├── covbind_results.csv
└── covbind_results.xlsx
```

## Requirements

The supplied `requirements.txt` contains:

``` text
requests>=2.31
pandas>=2.0
openpyxl>=3.1
python-dotenv>=1.0
```

The Python script also imports the `openai` package, but it is not
currently listed in the supplied `requirements.txt`.

Install the supplied dependencies:

``` bash
python -m pip install -r requirements.txt
```

Then install the additional package required by the script:

``` bash
python -m pip install openai
```

Alternatively, add this line to `requirements.txt`:

``` text
openai
```

and install everything with:

``` bash
python -m pip install -r requirements.txt
```

Python 3.9 or newer is recommended.

## Running the pipeline

The current script does not define command-line arguments. Its settings
are constants near the top of `covbind_pipeline_openai.py`.

Run it with:

``` bash
python covbind_pipeline_openai.py
```

The default output files are:

``` text
covbind_results.csv
covbind_results.xlsx
```

## Configuration

### Search query

`SEARCH_QUERY` defines the Europe PMC literature query.

The supplied query covers alpha- and betacoronaviruses and includes
terms such as SARS-CoV-2, SARS-CoV, MERS-CoV, HCoV-229E, HCoV-NL63,
HCoV-OC43, and HCoV-HKU1.

It combines those virus terms with protein-binding concepts such as
protein binding, receptor binding, protein-protein interaction, binding
domain, spike protein, and receptor-binding domain.

Edit `SEARCH_QUERY` in the Python file to change the literature scope.

### Maximum search results

``` python
MAX_RESULTS = 200
```

This limits the number of candidate Europe PMC records retrieved.

### Journal metric threshold

``` python
MIN_IMPACT_FACTOR = 4.0
```

The journal metric reader accepts a metric column named either `sjr` or
`if`. The configured threshold is applied to whichever supported metric
is present.

### Journal metric file

``` python
JOURNAL_METRIC_CSV = Path("scimagojr.csv")
```

The file is optional. If it is missing, the pipeline prints a warning
and continues without journal filtering.

### Output files

``` python
OUTPUT_XLSX = Path("covbind_results.xlsx")
OUTPUT_CSV = Path("covbind_results.csv")
```

Change these values to use different output names or locations.

### Model

``` python
OPENAI_MODEL = "o4-mini"
```

Change this value if a different model is required by the configured API
service.

## API configuration

The supplied script creates an OpenAI client using the custom API
endpoint:

``` text
https://apps.inside.anl.gov/argoapi/v1
```

The current source assigns a value directly to `OPENAI_API_KEY`.

For shared or version-controlled code, credentials should generally be
read from an environment variable instead of being stored in source
code. For example:

``` python
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "").strip()
```

Then set the environment variable before running the program.

macOS/Linux:

``` bash
export OPENAI_API_KEY="your-key-here"
python covbind_pipeline_openai.py
```

Windows PowerShell:

``` powershell
$env:OPENAI_API_KEY="your-key-here"
python covbind_pipeline_openai.py
```

Because the project includes `python-dotenv`, a `.env` file can also be
used after configuring the script to read the environment variable:

``` text
OPENAI_API_KEY=your-key-here
```

Do not commit `.env` files containing credentials to source control.

## Stage 1: Europe PMC search

The script queries the Europe PMC REST search API and requests core
metadata in JSON format.

For each paper it stores:

-   title
-   journal
-   publication year
-   DOI
-   abstract
-   Europe PMC source ID

Results are retrieved in pages of up to 100 records until `MAX_RESULTS`
is reached or no additional records are returned.

Internet access is required for this stage.

## Stage 2: Journal metric filtering

The pipeline attempts to read:

``` text
scimagojr.csv
```

The file may be comma-separated or semicolon-separated.

### Required columns

The journal-name column must be named either:

``` text
title
```

or:

``` text
journal
```

The metric column must be named either:

``` text
sjr
```

or:

``` text
if
```

Column-name detection is case-insensitive.

Example:

``` csv
Title,SJR
Journal of Example Virology,5.4
Another Journal,3.2
```

With the default threshold of `4.0`, the first example passes the
threshold and the second does not.

### Journal matching behavior

Journal names are converted to lowercase and then matched exactly.

Small differences in journal naming may therefore prevent a match. For
example:

``` text
Journal of Example Virology
```

and:

``` text
J. Example Virology
```

will not automatically be normalized to the same journal.

If the metric file does not exist, or does not contain the expected
columns, the pipeline continues without journal filtering.

## Stage 3: Structured extraction

For each paper containing an abstract, the model-based extractor
attempts to return four fields:

  -----------------------------------------------------------------------
  Field                               Description
  ----------------------------------- -----------------------------------
  `protein`                           Specific coronavirus protein
                                      discussed.

  `region`                            Region, domain, or residue range
                                      when stated.

  `binds_to`                          Host receptor or other binding
                                      partner.

  `importance`                        One or two sentences explaining the
                                      importance of the interaction based
                                      on the abstract.
  -----------------------------------------------------------------------

The extraction prompt instructs the model to return a JSON object and to
use `NA` when information is unavailable.

If the abstract is not about coronavirus protein binding, all four
extracted fields are expected to be `NA`.

## Heuristic fallback

If model extraction fails, the script uses a simple fallback.

The predefined protein vocabulary includes:

``` text
spike
nucleocapsid
ORF9b
ORF3a
envelope protein
membrane protein
nsp1
nsp3
nsp5
nsp12
nsp13
```

The predefined receptor/interaction vocabulary includes:

``` text
ACE2
DPP4
CD26
aminopeptidase N
APN
CD13
sialic acid
TOM70
heparan sulfate
TMPRSS2
```

The fallback also recognizes residue ranges such as:

``` text
residues 319-541
```

This heuristic is intentionally limited and is not equivalent to
model-based extraction.

## Paper filtering

A paper can be excluded from the final output when:

-   journal metrics are loaded and the journal does not meet the
    configured threshold;
-   the paper has no abstract; or
-   extraction produces both `protein = NA` and `binds_to = NA`.

The script prints summary counts at the end of the run.

## Output columns

Both output files contain the following columns:

  -----------------------------------------------------------------------
  Column                              Description
  ----------------------------------- -----------------------------------
  `Protein`                           Extracted viral protein name(s).

  `Region/Domain`                     Extracted protein region, domain,
                                      or residue range.

  `Binds To`                          Extracted receptor or interaction
                                      partner.

  `Importance`                        Abstract-based explanation of why
                                      the interaction matters.

  `Reference`                         Title, publication year, journal,
                                      and DOI when available.

  `Journal`                           Journal name from Europe PMC.

  `Year`                              Publication year.

  `DOI`                               Digital Object Identifier when
                                      available.

  `Journal Metric (SJR or IF)`        Matched journal metric when journal
                                      filtering is active.
  -----------------------------------------------------------------------

## Complete example

Install dependencies:

``` bash
python -m pip install -r requirements.txt
python -m pip install openai
```

Optionally place the journal metric file in the project directory:

``` text
scimagojr.csv
```

Then run:

``` bash
python covbind_pipeline_openai.py
```

A typical run prints messages such as:

``` text
[INFO] Searching Europe PMC (up to 200 results)...
[INFO] Retrieved ... candidate papers.
[INFO] Loaded metrics for ... journals from scimagojr.csv
[INFO] (1/...) Extracting: ...
[INFO] Done. Kept ... rows | dropped ... for journal metric | dropped ... for missing abstract.
[INFO] Wrote covbind_results.csv and covbind_results.xlsx
```

## Running without journal metrics

`scimagojr.csv` is optional.

If it is missing, the script reports:

``` text
[WARN] Journal metric file 'scimagojr.csv' not found.
       Continuing WITHOUT journal filtering -- all papers will be kept.
```

The literature search and extraction stages continue normally.

## Troubleshooting

### `ModuleNotFoundError: No module named 'openai'`

The Python file imports `OpenAI`, but `openai` is not included in the
supplied `requirements.txt`.

Install it with:

``` bash
python -m pip install openai
```

### Other missing Python packages

Run:

``` bash
python -m pip install -r requirements.txt
```

### Excel output fails

Excel export requires `openpyxl`:

``` bash
python -m pip install openpyxl
```

### Europe PMC requests fail

Check internet connectivity, proxy/firewall settings, and Europe PMC
availability. The script uses a 30-second HTTP timeout.

### Model extraction fails

The script catches model-extraction errors and falls back to heuristic
extraction. Look for:

``` text
[WARN] OpenAI extraction failed (...); falling back to heuristic.
```

### No papers pass the journal filter

Possible causes include:

-   the threshold is too high;
-   the journal is missing from the metric CSV;
-   Europe PMC and the metric file use different journal-name variants;
    or
-   the metric file does not contain the expected columns.

### Output contains few records

The pipeline applies multiple filters. Papers can be excluded because of
journal metrics, missing abstracts, or extraction that does not identify
both a protein or binding partner.

## Reproducibility

Pipeline results can change as the literature database, journal metric
data, search query, extraction model, and API behavior change.

For reproducible analyses, retain:

-   the exact Python script;
-   `requirements.txt`;
-   the journal metric CSV used for the run;
-   generated CSV and Excel files;
-   model/configuration information; and
-   the date of the run.

## Interpretation and limitations

The output is an automated literature-screening and extraction aid.

Important limitations include:

-   extraction is based on abstracts rather than necessarily full text;
-   journal matching uses exact lowercased journal names;
-   the heuristic fallback has a fixed vocabulary;
-   journal metrics do not establish the validity of an individual paper
    or extracted claim; and
-   extracted interactions should be checked against the cited
    publication before use in downstream scientific conclusions.

## License

No license is specified in the supplied project files. Add an
appropriate license if the project will be distributed.

## Notes
- Authored by Omar Loay, co-authored by Claude and Gemini based on the script "viral_protein_pipeline.py" by Milana Djonovic. Modifications made by Sachin Kumar (include base url for Argo server and provide API_KEY)
