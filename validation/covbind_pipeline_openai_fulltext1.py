#!/usr/bin/env python3
"""
covbind_pipeline.py (OpenAI version, Semantic Scholar search + full-text XML)

A small, self-contained literature pipeline for finding papers about
alpha- and betacoronavirus PROTEIN BINDING, restricted to journals above a
minimum impact threshold, using OpenAI for structured extraction.

CHANGE FROM PREVIOUS VERSION:
  - Search/metadata now comes from the Semantic Scholar Graph API instead
    of Europe PMC.
  - Full text is now fetched when available. Semantic Scholar does NOT
    host full-text XML itself, so for full text we take the PMC ID that
    Semantic Scholar returns (externalIds.PubMedCentral) and fetch the
    actual XML from Europe PMC's fullTextXML endpoint (this only works
    for papers in the PMC Open Access subset). If that fails or no PMC ID
    exists, the pipeline falls back to the abstract Semantic Scholar gave
    us -- no paper is silently dropped just because full text isn't
    available.

Outputs: covbind_results.xlsx and covbind_results.csv

-----------------------------------------------------------------------
SETUP NOTES
-----------------------------------------------------------------------
- Semantic Scholar's public API works without a key, but with low rate
  limits (roughly 1 request / few seconds, shared across all anonymous
  users). If you hit 429s, get a free API key at
  https://www.semanticscholar.org/product/api and set:
      export S2_API_KEY="your-key-here"
  The script will send it automatically if present.
- No new dependency for the XML step -- it uses Python's built-in
  xml.etree.ElementTree.
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import requests
import pandas as pd

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from openai import OpenAI

# =========================================================================
# CONFIG
# =========================================================================

# Semantic Scholar's search syntax is simpler than Europe PMC's -- it does
# a relevance-ranked text search rather than a field-query boolean, so we
# pass it as a plain query string (S2 tokenizes/matches internally). You
# can still bias it with quoted phrases.
SEARCH_QUERY = (
    'alphacoronavirus betacoronavirus SARS-CoV-2 SARS-CoV MERS-CoV '
    'HCoV-229E HCoV-NL63 HCoV-OC43 HCoV-HKU1 protein binding receptor '
    'binding domain spike protein protein-protein interaction'
)

MAX_RESULTS = 200
MIN_IMPACT_FACTOR = 4.0

JOURNAL_METRIC_CSV = Path("scimagojr.csv")
OUTPUT_XLSX = Path("covbind_results.xlsx")
OUTPUT_CSV = Path("covbind_results.csv")

SEMANTIC_SCHOLAR_SEARCH_URL = "https://api.semanticscholar.org/graph/v1/paper/search"
S2_API_KEY = "put your key here"
S2_FIELDS = "title,abstract,year,venue,journal,externalIds"

EUROPEPMC_FULLTEXT_TMPL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/{ext_id}/fullTextXML"
)

REQUEST_DELAY_SEC = 1.1  # Semantic Scholar's anonymous rate limit is tight
S2_MAX_RETRIES = 6       # give up after this many consecutive 429s on one request
S2_BACKOFF_BASE_SEC = 15 # first backoff wait; doubles each retry (capped below)
S2_BACKOFF_CAP_SEC = 180

USE_FULL_TEXT = True
FULLTEXT_MAX_CHARS = 15000
SECTION_KEEP_KEYWORDS = ("method", "result", "discussion")

OPENAI_MODEL = "GPT-5.6 Terra"
OPENAI_API_KEY = "ac.oloay"  # os.environ.get("OPENAI_API_KEY", "").strip()

client = OpenAI(api_key=OPENAI_API_KEY, base_url="https://apps.inside.anl.gov/argoapi/v1") if OPENAI_API_KEY else None


# =========================================================================
# DATA MODEL
# =========================================================================

@dataclass
class Paper:
    title: str
    journal: str
    year: str
    doi: str
    abstract: str
    source_id: str
    pmcid: str = ""  # e.g. "PMC1234567" -- only set when Semantic Scholar has it

@dataclass
class ExtractedRow:
    protein: str = "NA"
    region: str = "NA"
    binds_to: str = "NA"
    importance: str = "NA"
    reference: str = ""
    journal: str = ""
    year: str = ""
    doi: str = ""
    metric_value: Optional[float] = None
    text_source: str = ""  # 'europepmc_fulltext' | 'abstract'


# =========================================================================
# STAGE 1 -- SEARCH (Semantic Scholar Graph API)
# =========================================================================

def search_semanticscholar(query: str, max_results: int) -> list[Paper]:
    """
    Query the Semantic Scholar /paper/search endpoint and page through
    results via its 'next' offset cursor.
    Docs: https://api.semanticscholar.org/api-docs/graph#tag/Paper-Data/operation/get_graph_paper_relevance_search

    Rate-limit handling: on a 429, honors the Retry-After header if the
    API sends one, otherwise backs off exponentially starting at
    S2_BACKOFF_BASE_SEC and capping at S2_BACKOFF_CAP_SEC. Gives up after
    S2_MAX_RETRIES consecutive 429s on the SAME request (rather than
    retrying forever) and raises a clear error telling you to get an
    S2_API_KEY, since that's the actual fix for persistent 429s.
    """
    papers: list[Paper] = []
    offset = 0
    headers = {"x-api-key": S2_API_KEY} if S2_API_KEY else {}
    if not S2_API_KEY:
        print(
            "[WARN] No S2_API_KEY set. Semantic Scholar's anonymous tier is "
            "shared across all unauthenticated users and gets throttled hard "
            "-- if you see repeated 429s below, get a free key at "
            "https://www.semanticscholar.org/product/api and set "
            "S2_API_KEY, then re-run."
        )

    while len(papers) < max_results:
        params = {
            "query": query,
            "fields": S2_FIELDS,
            "offset": offset,
            "limit": min(100, max_results - len(papers)),  # S2 caps at 100/page
        }

        retries = 0
        while True:
            resp = requests.get(SEMANTIC_SCHOLAR_SEARCH_URL, params=params, headers=headers, timeout=30)
            if resp.status_code != 429:
                break

            retries += 1
            if retries > S2_MAX_RETRIES:
                raise RuntimeError(
                    f"Semantic Scholar kept returning 429 after {S2_MAX_RETRIES} "
                    f"retries. This almost always means the anonymous rate limit "
                    f"is exhausted right now. Get a free API key at "
                    f"https://www.semanticscholar.org/product/api, "
                    f"`export S2_API_KEY=...`, and re-run -- or just wait a few "
                    f"minutes and try again without a key."
                )

            retry_after = resp.headers.get("Retry-After")
            if retry_after and retry_after.strip().isdigit():
                wait_sec = int(retry_after.strip())
            else:
                wait_sec = min(S2_BACKOFF_BASE_SEC * (2 ** (retries - 1)), S2_BACKOFF_CAP_SEC)

            print(f"[WARN] Semantic Scholar 429 (attempt {retries}/{S2_MAX_RETRIES}); waiting {wait_sec}s...")
            time.sleep(wait_sec)

        resp.raise_for_status()  # raise on any non-429 error status (e.g. 5xx, 400)
        data = resp.json()

        results = data.get("data", [])
        if not results:
            break

        for r in results:
            ext = r.get("externalIds") or {}
            doi = (ext.get("DOI") or "").strip()

            raw_pmc = (ext.get("PubMedCentral") or "").strip()
            pmcid = raw_pmc if raw_pmc.upper().startswith("PMC") else (f"PMC{raw_pmc}" if raw_pmc else "")

            journal_obj = r.get("journal") or {}
            journal_name = (journal_obj.get("name") or "").strip() if isinstance(journal_obj, dict) else ""
            if not journal_name:
                journal_name = (r.get("venue") or "").strip()  # fallback, e.g. for preprints/conferences

            papers.append(Paper(
                title=(r.get("title") or "").strip(),
                journal=journal_name,
                year=str(r.get("year") or "").strip(),
                doi=doi,
                abstract=(r.get("abstract") or "").strip(),
                source_id=str(r.get("paperId") or ""),
                pmcid=pmcid,
            ))

        next_offset = data.get("next")
        if next_offset is None:
            break
        offset = next_offset
        time.sleep(REQUEST_DELAY_SEC)

    return papers[:max_results]


# =========================================================================
# STAGE 1b -- FULL TEXT RETRIEVAL (Europe PMC XML, keyed off the PMC ID
# that Semantic Scholar gave us)
# =========================================================================

# Diagnostic counters so you can see WHY papers fall back to abstract,
# instead of just seeing that they did. Printed at the end of run().
FULLTEXT_STATS = {
    "no_pmcid": 0,        # Semantic Scholar gave no PMC ID at all
    "http_error": 0,      # PMC ID exists, but Europe PMC returned non-200
    "parse_error": 0,     # got a 200, but the XML didn't parse
    "empty_sections": 0,  # parsed fine, but no usable sec text found
    "success": 0,
}


def fetch_fulltext_epmc(paper: Paper, verbose: bool = False) -> Optional[str]:
    """
    Fetch full-text XML from Europe PMC using the PMC ID from Semantic
    Scholar's externalIds. Only works for papers in the PMC Open Access
    subset -- returns None for anything else (paywalled, no PMC ID, etc.),
    and the caller falls back to the abstract in that case.

    Updates FULLTEXT_STATS so the run summary can tell you which failure
    mode is actually happening, instead of everything looking identical
    as "fell back to abstract".
    """
    if not paper.pmcid:
        FULLTEXT_STATS["no_pmcid"] += 1
        return None

    url = EUROPEPMC_FULLTEXT_TMPL.format(ext_id=paper.pmcid)
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200 or not resp.content:
            FULLTEXT_STATS["http_error"] += 1
            if verbose:
                print(f"    [fulltext] {paper.pmcid}: HTTP {resp.status_code}")
            return None
        root = ET.fromstring(resp.content)
    except requests.RequestException as e:
        FULLTEXT_STATS["http_error"] += 1
        if verbose:
            print(f"    [fulltext] {paper.pmcid}: request failed ({e})")
        return None
    except ET.ParseError as e:
        FULLTEXT_STATS["parse_error"] += 1
        if verbose:
            print(f"    [fulltext] {paper.pmcid}: XML parse failed ({e})")
        return None

    chunks = []
    for sec in root.iter("sec"):
        title_el = sec.find("title")
        title_text = (title_el.text or "").lower() if title_el is not None else ""
        if any(k in title_text for k in SECTION_KEEP_KEYWORDS) or not chunks:
            text = " ".join(t.strip() for t in sec.itertext() if t and t.strip())
            if text:
                chunks.append(text)

    full_text = "\n\n".join(chunks).strip()
    if not full_text:
        FULLTEXT_STATS["empty_sections"] += 1
        if verbose:
            print(f"    [fulltext] {paper.pmcid}: XML parsed but no <sec> text found")
        return None

    FULLTEXT_STATS["success"] += 1
    return full_text


def get_source_text(paper: Paper, verbose: bool = False) -> tuple[str, str]:
    """Returns (text, method) -- 'europepmc_fulltext' or 'abstract'."""
    if USE_FULL_TEXT:
        text = fetch_fulltext_epmc(paper, verbose=verbose)
        if text:
            return text[:FULLTEXT_MAX_CHARS], "europepmc_fulltext"
    return paper.abstract, "abstract"


# =========================================================================
# STAGE 2 -- JOURNAL QUALITY GATE
# =========================================================================

def load_journal_metrics(csv_path: Path) -> dict[str, float]:
    if not csv_path.exists():
        print(
            f"[WARN] Journal metric file '{csv_path}' not found.\n"
            f"       Continuing WITHOUT journal filtering -- all papers will be kept.\n"
        )
        return {}

    lookup: dict[str, float] = {}

    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)
        f.seek(0)
        delimiter = ";" if sample.count(";") > sample.count(",") else ","
        reader = csv.DictReader(f, delimiter=delimiter)

        fieldnames = [fn.strip().lower() for fn in (reader.fieldnames or [])]
        title_col = None
        metric_col = None
        for candidate in ("title", "journal"):
            if candidate in fieldnames:
                title_col = reader.fieldnames[fieldnames.index(candidate)]
                break
        for candidate in ("sjr", "if"):
            if candidate in fieldnames:
                metric_col = reader.fieldnames[fieldnames.index(candidate)]
                break

        if not title_col or not metric_col:
            print(f"[WARN] Missing columns in {csv_path}. Skipping journal filter.")
            return {}

        for row in reader:
            name = (row.get(title_col) or "").strip().lower()
            raw_val = (row.get(metric_col) or "").strip().replace(",", ".")
            if not name or not raw_val:
                continue
            try:
                lookup[name] = float(raw_val)
            except ValueError:
                continue

    print(f"[INFO] Loaded metrics for {len(lookup)} journals from {csv_path}")
    return lookup


def journal_passes(journal_name: str, metrics: dict[str, float], min_value: float) -> Optional[float]:
    if not metrics:
        return None
    val = metrics.get(journal_name.strip().lower())
    if val is not None and val >= min_value:
        return val
    return None


# =========================================================================
# STAGE 3 -- EXTRACTION (OpenAI with Heuristic Fallback)
# =========================================================================

EXTRACTION_PROMPT = """You are extracting structured facts from a virology \
paper excerpt (this may be an abstract, or full Methods/Results/Discussion \
text) about coronavirus PROTEIN BINDING (e.g. viral protein binding host \
receptor, or protein-protein interactions).

Return ONLY a JSON object with these exact keys:
  "protein":     the specific protein discussed (e.g. "Spike (S) protein", "ORF9b"). If multiple, comma-separate.
  "region":      the specific region/domain/residue range, if any (e.g. "RBD, residues 319-541"). Use "NA" if not specified.
  "binds_to":    what the protein binds (e.g. "human ACE2", "DPP4"). Use "NA" if unclear.
  "importance":  1-2 sentences on why this interaction matters based ONLY on the text. Use "NA" if not stated.

If the text is not about coronavirus protein binding, set all four fields to "NA".

Text:
\"\"\"{text}\"\"\"
"""


def extract_with_openai(text: str) -> dict:
    """Call OpenAI API for structured extraction using JSON mode."""
    response = client.chat.completions.create(
        model=OPENAI_MODEL,
        response_format={"type": "json_object"},
        messages=[
            {"role": "system", "content": "You extract structured data from literature text into JSON."},
            {"role": "user", "content": EXTRACTION_PROMPT.format(text=text)}
        ]
    )
    return json.loads(response.choices[0].message.content)


KNOWN_PROTEINS = [
    "spike", "nucleocapsid", "ORF9b", "ORF3a", "envelope protein",
    "membrane protein", "nsp1", "nsp3", "nsp5", "nsp12", "nsp13",
]
KNOWN_RECEPTORS = [
    "ACE2", "DPP4", "CD26", "aminopeptidase N", "APN", "CD13",
    "sialic acid", "TOM70", "heparan sulfate", "TMPRSS2",
]

def extract_heuristic(text: str) -> dict:
    text_lower = text.lower()
    proteins = [p for p in KNOWN_PROTEINS if p.lower() in text_lower]
    receptors = [r for r in KNOWN_RECEPTORS if r.lower() in text_lower]

    region = "NA"
    m = re.search(r"residues?\s+\d+\s*[-–]\s*\d+", text, flags=re.IGNORECASE)
    if m:
        region = m.group(0)

    return {
        "protein": ", ".join(proteins) if proteins else "NA",
        "region": region,
        "binds_to": ", ".join(receptors) if receptors else "NA",
        "importance": "NA (OPENAI_API_KEY not set - fallback heuristic used)",
    }


def extract_fields(text: str) -> dict:
    if not text:
        return {"protein": "NA", "region": "NA", "binds_to": "NA", "importance": "NA"}
    if client:
        try:
            return extract_with_openai(text)
        except Exception as e:
            print(f"[WARN] OpenAI extraction failed ({e}); falling back to heuristic.")
    return extract_heuristic(text)


# =========================================================================
# STAGE 4 -- MAIN PIPELINE
# =========================================================================

def build_reference(paper: Paper) -> str:
    year = paper.year or "n.d."
    journal = paper.journal or "Unknown journal"
    doi_part = f" https://doi.org/{paper.doi}" if paper.doi else ""
    return f"{paper.title} ({year}). {journal}.{doi_part}"


def run() -> None:
    print(f"[INFO] Searching Semantic Scholar (up to {MAX_RESULTS} results)...")
    papers = search_semanticscholar(SEARCH_QUERY, MAX_RESULTS)
    print(f"[INFO] Retrieved {len(papers)} candidate papers.")

    metrics = load_journal_metrics(JOURNAL_METRIC_CSV)

    rows: list[ExtractedRow] = []
    kept, dropped_no_journal_match, dropped_no_text = 0, 0, 0

    for i, paper in enumerate(papers, start=1):
        metric_value = journal_passes(paper.journal, metrics, MIN_IMPACT_FACTOR)

        if metrics and metric_value is None:
            dropped_no_journal_match += 1
            continue

        source_text, text_source = get_source_text(paper, verbose=True)
        if not source_text:
            dropped_no_text += 1
            continue

        print(f"[INFO] ({i}/{len(papers)}) Extracting [{text_source}]: {paper.title[:70]}...")
        fields = extract_fields(source_text)

        if fields.get("protein", "NA") == "NA" and fields.get("binds_to", "NA") == "NA":
            continue

        rows.append(ExtractedRow(
            protein=fields.get("protein", "NA"),
            region=fields.get("region", "NA"),
            binds_to=fields.get("binds_to", "NA"),
            importance=fields.get("importance", "NA"),
            reference=build_reference(paper),
            journal=paper.journal,
            year=paper.year,
            doi=paper.doi,
            metric_value=metric_value,
            text_source=text_source,
        ))
        kept += 1

        if client:
            time.sleep(0.2)

    print(
        f"[INFO] Done. Kept {kept} rows | "
        f"dropped {dropped_no_journal_match} for journal metric | "
        f"dropped {dropped_no_text} for no usable text (full text or abstract)."
    )
    print(
        "[INFO] Full-text retrieval breakdown: "
        f"no PMC ID={FULLTEXT_STATS['no_pmcid']}, "
        f"HTTP/request error={FULLTEXT_STATS['http_error']}, "
        f"XML parse error={FULLTEXT_STATS['parse_error']}, "
        f"parsed but empty={FULLTEXT_STATS['empty_sections']}, "
        f"success={FULLTEXT_STATS['success']}"
    )

    df = pd.DataFrame([asdict(r) for r in rows])
    df = df.rename(columns={
        "protein": "Protein",
        "region": "Region/Domain",
        "binds_to": "Binds To",
        "importance": "Importance",
        "reference": "Reference",
        "journal": "Journal",
        "year": "Year",
        "doi": "DOI",
        "metric_value": "Journal Metric (SJR or IF)",
        "text_source": "Extracted From",
    })

    df.to_csv(OUTPUT_CSV, index=False)
    df.to_excel(OUTPUT_XLSX, index=False)
    print(f"[INFO] Wrote {OUTPUT_CSV} and {OUTPUT_XLSX}")


if __name__ == "__main__":
    run()
