#!/usr/bin/env python3
"""
covbind_pipeline.py (OpenAI Version)

A small, self-contained literature pipeline for finding papers about
alpha- and betacoronavirus PROTEIN BINDING, restricted to journals above a
minimum impact threshold, using OpenAI for structured extraction.

Outputs: covbind_results.xlsx and covbind_results.csv
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
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

SEARCH_QUERY = (
    '(alphacoronavirus OR betacoronavirus OR "alpha coronavirus" OR '
    '"beta coronavirus" OR SARS-CoV-2 OR SARS-CoV OR MERS-CoV OR '
    '"HCoV-229E" OR "HCoV-NL63" OR "HCoV-OC43" OR "HCoV-HKU1") '
    'AND ("protein binding" OR "receptor binding" OR "protein-protein '
    'interaction" OR "binding domain" OR "spike protein" OR "receptor '
    'binding domain")'
)

MAX_RESULTS = 200
MIN_IMPACT_FACTOR = 4.0

JOURNAL_METRIC_CSV = Path("scimagojr.csv")
OUTPUT_XLSX = Path("covbind_results.xlsx")
OUTPUT_CSV = Path("covbind_results.csv")

EUROPEPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
REQUEST_DELAY_SEC = 0.34

OPENAI_MODEL = "o4-mini"
OPENAI_API_KEY = "ac.skumar"#os.environ.get("OPENAI_API_KEY", "").strip()

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


# =========================================================================
# STAGE 1 -- SEARCH (Europe PMC)
# =========================================================================

def search_europepmc(query: str, max_results: int) -> list[Paper]:
    papers: list[Paper] = []
    cursor_mark = "*"
    page_size = 100

    while len(papers) < max_results:
        params = {
            "query": query,
            "format": "json",
            "pageSize": min(page_size, max_results - len(papers)),
            "cursorMark": cursor_mark,
            "resultType": "core",
        }
        resp = requests.get(EUROPEPMC_BASE, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        results = data.get("resultList", {}).get("result", [])
        if not results:
            break

        for r in results:
            papers.append(Paper(
                title=r.get("title", "").strip(),
                journal=r.get("journalTitle", "").strip(),
                year=str(r.get("pubYear", "")).strip(),
                doi=r.get("doi", "").strip(),
                abstract=r.get("abstractText", "").strip(),
                source_id=r.get("id", ""),
            ))

        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor_mark:
            break
        cursor_mark = next_cursor
        time.sleep(REQUEST_DELAY_SEC)

    return papers[:max_results]


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
abstract about coronavirus PROTEIN BINDING (e.g. viral protein binding host \
receptor, or protein-protein interactions).

Return ONLY a JSON object with these exact keys:
  "protein":     the specific protein discussed (e.g. "Spike (S) protein", "ORF9b"). If multiple, comma-separate.
  "region":      the specific region/domain/residue range, if any (e.g. "RBD, residues 319-541"). Use "NA" if not specified.
  "binds_to":    what the protein binds (e.g. "human ACE2", "DPP4"). Use "NA" if unclear.
  "importance":  1-2 sentences on why this interaction matters based ONLY on the abstract. Use "NA" if not stated.

If the abstract is not about coronavirus protein binding, set all four fields to "NA".

Abstract:
\"\"\"{abstract}\"\"\"
"""


def extract_with_openai(abstract: str) -> dict:
    """Call OpenAI API for structured extraction using JSON mode."""
    response = client.chat.completions.create(
        model=OPENAI_MODEL,
        response_format={"type": "json_object"},
        messages=[
            {"role": "system", "content": "You extract structured data from literature abstracts into JSON."},
            {"role": "user", "content": EXTRACTION_PROMPT.format(abstract=abstract)}
        ],
        temperature=0.0
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

def extract_heuristic(abstract: str) -> dict:
    text_lower = abstract.lower()
    proteins = [p for p in KNOWN_PROTEINS if p.lower() in text_lower]
    receptors = [r for r in KNOWN_RECEPTORS if r.lower() in text_lower]

    region = "NA"
    m = re.search(r"residues?\s+\d+\s*[-–]\s*\d+", abstract, flags=re.IGNORECASE)
    if m:
        region = m.group(0)

    return {
        "protein": ", ".join(proteins) if proteins else "NA",
        "region": region,
        "binds_to": ", ".join(receptors) if receptors else "NA",
        "importance": "NA (OPENAI_API_KEY not set - fallback heuristic used)",
    }


def extract_fields(abstract: str) -> dict:
    if not abstract:
        return {"protein": "NA", "region": "NA", "binds_to": "NA", "importance": "NA"}
    if client:
        try:
            return extract_with_openai(abstract)
        except Exception as e:
            print(f"[WARN] OpenAI extraction failed ({e}); falling back to heuristic.")
    return extract_heuristic(abstract)


# =========================================================================
# STAGE 4 -- MAIN PIPELINE
# =========================================================================

def build_reference(paper: Paper) -> str:
    year = paper.year or "n.d."
    journal = paper.journal or "Unknown journal"
    doi_part = f" https://doi.org/{paper.doi}" if paper.doi else ""
    return f"{paper.title} ({year}). {journal}.{doi_part}"


def run() -> None:
    print(f"[INFO] Searching Europe PMC (up to {MAX_RESULTS} results)...")
    papers = search_europepmc(SEARCH_QUERY, MAX_RESULTS)
    print(f"[INFO] Retrieved {len(papers)} candidate papers.")

    metrics = load_journal_metrics(JOURNAL_METRIC_CSV)

    rows: list[ExtractedRow] = []
    kept, dropped_no_journal_match, dropped_no_abstract = 0, 0, 0

    for i, paper in enumerate(papers, start=1):
        metric_value = journal_passes(paper.journal, metrics, MIN_IMPACT_FACTOR)

        if metrics and metric_value is None:
            dropped_no_journal_match += 1
            continue

        if not paper.abstract:
            dropped_no_abstract += 1
            continue

        print(f"[INFO] ({i}/{len(papers)}) Extracting: {paper.title[:80]}...")
        fields = extract_fields(paper.abstract)

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
        ))
        kept += 1

        if client:
            time.sleep(0.2)

    print(
        f"[INFO] Done. Kept {kept} rows | "
        f"dropped {dropped_no_journal_match} for journal metric | "
        f"dropped {dropped_no_abstract} for missing abstract."
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
    })

    df.to_csv(OUTPUT_CSV, index=False)
    df.to_excel(OUTPUT_XLSX, index=False)
    print(f"[INFO] Wrote {OUTPUT_CSV} and {OUTPUT_XLSX}")


if __name__ == "__main__":
    run()