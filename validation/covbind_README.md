# CovBind Literature Pipeline — Semantic Scholar + Full-Text Extraction

`covbind_pipeline_openai_fulltext1.py` searches for alpha- and betacoronavirus protein-binding literature, optionally filters papers by journal metric, retrieves full text when available, extracts structured binding information with an OpenAI-compatible model, and exports the results to CSV and Excel.

## What changed in this version

This version uses the **Semantic Scholar Graph API** for search and metadata instead of Europe PMC search. When Semantic Scholar supplies a PubMed Central ID, the script attempts to retrieve full-text XML from Europe PMC. If full text is unavailable, it falls back to the Semantic Scholar abstract.

Default outputs:

```text
covbind_results.csv
covbind_results.xlsx
```

## Workflow

```text
Semantic Scholar search
        |
        v
candidate papers
        |
        v
optional journal metric filter
        |
        v
PMC ID available?
   |             |
  yes            no
   |             |
   v             |
Europe PMC       |
fullTextXML      |
   |             |
   +------> usable full text?
              |       |
             yes      no
              |       |
              v       v
          full text  abstract
              \       /
               \     /
                v   v
          structured extraction
                |
                v
      covbind_results.csv/.xlsx
```

## Requirements

Python 3.9+ is recommended.

Install:

```bash
python -m pip install requests pandas openpyxl python-dotenv openai
```

Suggested `requirements.txt`:

```text
requests>=2.31
pandas>=2.0
openpyxl>=3.1
python-dotenv>=1.0
openai
```

The full-text XML step uses Python's built-in `xml.etree.ElementTree`, so no additional XML package is required.

## Running the script

The current script has no command-line arguments. Configure the constants near the top of the file and run:

```bash
python covbind_pipeline_openai_fulltext1.py
```

## Main configuration

```python
MAX_RESULTS = 200
MIN_IMPACT_FACTOR = 4.0

JOURNAL_METRIC_CSV = Path("scimagojr.csv")
OUTPUT_XLSX = Path("covbind_results.xlsx")
OUTPUT_CSV = Path("covbind_results.csv")

USE_FULL_TEXT = True
FULLTEXT_MAX_CHARS = 15000
SECTION_KEEP_KEYWORDS = ("method", "result", "discussion")

OPENAI_MODEL = "GPT-5.6 Terra"
```

`SEARCH_QUERY` controls the Semantic Scholar search terms.

## Semantic Scholar API

The search stage requests paper metadata including title, abstract, year, venue/journal, and external identifiers.

The script supports a Semantic Scholar API key. For shared code, prefer reading it from the environment:

```python
S2_API_KEY = os.environ.get("S2_API_KEY", "").strip()
```

Then:

```bash
export S2_API_KEY="your-key-here"
```

Without a key, anonymous Semantic Scholar access may be rate-limited.

### Rate-limit handling

The script handles HTTP 429 responses using `Retry-After` when supplied and otherwise exponential backoff. Relevant settings include:

```python
REQUEST_DELAY_SEC = 1.1
S2_MAX_RETRIES = 6
S2_BACKOFF_BASE_SEC = 15
S2_BACKOFF_CAP_SEC = 180
```

## Full-text retrieval

When `USE_FULL_TEXT = True`, the script uses the PMC ID returned by Semantic Scholar to request Europe PMC full-text XML.

If full text cannot be used because there is no PMC ID, an HTTP/request failure occurs, XML parsing fails, or no usable section text is found, the paper falls back to its abstract.

The script reports diagnostic counts for:

```text
no_pmcid
http_error
parse_error
empty_sections
success
```

Full-text text passed to extraction is limited by:

```python
FULLTEXT_MAX_CHARS = 15000
```

The XML processing preferentially retains sections whose titles include `method`, `result`, or `discussion`.

## Journal metric filtering

The optional journal metric file defaults to:

```text
scimagojr.csv
```

Supported journal-name columns:

```text
title
journal
```

Supported metric columns:

```text
sjr
if
```

Column matching is case-insensitive. Journal-name matching itself is lowercased exact matching.

If the metric file is missing or lacks usable columns, the script continues without journal filtering.

## Structured extraction

The model extracts:

| Field | Description |
| --- | --- |
| `protein` | Coronavirus protein discussed. |
| `region` | Region/domain/residue range when stated. |
| `binds_to` | Receptor or other binding partner. |
| `importance` | Why the interaction matters, based only on the supplied text. |

The source can be either a Europe PMC full-text excerpt or the Semantic Scholar abstract.

If the text is not about coronavirus protein binding, the extraction prompt requests `NA` for all four fields.

## Heuristic fallback

If model extraction fails, the script falls back to keyword/regex extraction.

Recognized protein terms include spike, nucleocapsid, ORF9b, ORF3a, envelope protein, membrane protein, nsp1, nsp3, nsp5, nsp12, and nsp13.

Recognized interaction terms include ACE2, DPP4, CD26, aminopeptidase N, APN, CD13, sialic acid, TOM70, heparan sulfate, and TMPRSS2.

The fallback can also identify residue ranges such as:

```text
residues 319-541
```

## Output columns

`covbind_results.csv` and `covbind_results.xlsx` contain:

```text
Protein
Region/Domain
Binds To
Importance
Reference
Journal
Year
DOI
Journal Metric (SJR or IF)
Extracted From
```

`Extracted From` identifies whether extraction used:

```text
europepmc_fulltext
```

or:

```text
abstract
```

## Example setup

```bash
python -m venv .venv
source .venv/bin/activate

python -m pip install requests pandas openpyxl python-dotenv openai

export S2_API_KEY="your-semantic-scholar-key"
export OPENAI_API_KEY="your-api-key"

python covbind_pipeline_openai_fulltext1.py
```

On Windows PowerShell, activation/environment-variable syntax will differ.

## Credential configuration

The supplied script contains API configuration in source. For shared or version-controlled use, change credentials to environment-based configuration, for example:

```python
S2_API_KEY = os.environ.get("S2_API_KEY", "").strip()
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "").strip()
```

A `.env` file may then contain:

```text
S2_API_KEY=your-semantic-scholar-key
OPENAI_API_KEY=your-api-key
```

Add `.env` to `.gitignore` and do not commit credentials.

## Troubleshooting

### Semantic Scholar returns 429

Configure `S2_API_KEY`, reduce request volume, or retry later. The script already implements bounded retry/backoff behavior.

### Full text frequently falls back to abstract

This is expected when papers lack a PMC ID or Europe PMC full-text XML is unavailable. Check the full-text diagnostic counters printed at the end.

### `scimagojr.csv` is ignored

Confirm that the file contains a `title` or `journal` column and an `sjr` or `if` column.

### `ModuleNotFoundError: openai`

```bash
python -m pip install openai
```

### Excel export fails

```bash
python -m pip install openpyxl
```

### Model extraction fails

The script prints a warning and uses the heuristic fallback.

## Reproducibility and limitations

For reproducible analyses, retain the exact script, dependency versions, journal metric file, generated outputs, search query, model configuration, and run date.

Important limitations include:

- Semantic Scholar search is relevance-ranked rather than an exhaustive fielded systematic search.
- Full text is available only for papers for which suitable Europe PMC XML can be retrieved.
- Full-text input is truncated before extraction.
- Journal-name matching is exact after lowercasing.
- The heuristic fallback uses a fixed vocabulary.
- Automated extracted claims should be checked against the cited publication before downstream scientific interpretation.
- Authored by Omar Loay, using Claude, GPT and Gemini.