#!/usr/bin/env python3
"""
Viral Protein Processing Pipeline
Retrieves, filters (QC), annotates, aligns, and identifies conserved regions
in viral protein sequences from NCBI.
"""

import pandas as pd
import os
import sys
import json
import re
import math
import argparse
import urllib.request
import urllib.parse
import time
from typing import List, Dict, Any, Tuple

# NCBI E-utilities Base URLs
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

def http_get(url: str, params: Dict[str, Any]) -> str:
    """Helper to perform HTTP GET with rate limiting and exponential backoff."""
    encoded_params = urllib.parse.urlencode(params)
    full_url = f"{url}?{encoded_params}"
    
    retries = 3
    backoff = 2
    for attempt in range(retries):
        try:
            req = urllib.request.Request(full_url, headers={'User-Agent': 'AntigravityViralPipeline/1.0'})
            with urllib.request.urlopen(req) as response:
                return response.read().decode('utf-8')
        except Exception as e:
            if attempt == retries - 1:
                raise RuntimeError(f"HTTP request failed for {url}: {e}")
            time.sleep(backoff)
            backoff *= 2
    return ""

# =====================================================================
# STEP 1: RETRIEVAL
# =====================================================================
def fetch_viral_proteins(taxid: str, limit: int, output_file: str, host_filter: bool = False):
    print(f"[*] Searching NCBI Entrez for viral proteins (TaxID: {taxid}, Limit: {limit}, HostFilter: {host_filter})...")
    
    # Base search term
    term = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP] AND protein[biomaterial/molecule type]"
    if host_filter:
        # Search for host keywords: mammals, human, bat, swine, bovine, camel, pangolin, etc.
        host_term = '(mammal[Host] OR pangolin[Host] OR human[Host] OR bat[Host] OR swine[Host] OR bovine[Host] OR camel[Host] OR "Manis javanica"[Organism])'
        term = f"{term} AND {host_term}"
        
    esearch_params = {
        "db": "protein",
        "term": term,
        "retmode": "json",
        "retmax": limit
    }
    
    response_text = http_get(ESEARCH_URL, esearch_params)
    search_data = json.loads(response_text)
    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    count = search_data.get("esearchresult", {}).get("count", "0")
    
    # Fallback if host_term filter yielded 0 results due to unindexed host fields in NCBI protein
    if not id_list and host_filter:
        print("[!] No records matched strict host index. Falling back to Entrez text query for mammals/pangolins...")
        term = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP] AND protein[biomaterial/molecule type] AND (mammal OR pangolin OR bat OR human OR swine OR bovine OR camel)"
        esearch_params["term"] = term
        response_text = http_get(ESEARCH_URL, esearch_params)
        search_data = json.loads(response_text)
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        count = search_data.get("esearchresult", {}).get("count", "0")

    print(f"[*] Found {count} matching records. Fetching top {len(id_list)} records...")
    
    if not id_list:
        print("[!] No records found.")
        records = []
    else:
        esummary_params = {
            "db": "protein",
            "id": ",".join(id_list),
            "retmode": "json"
        }
        summary_text = http_get(ESUMMARY_URL, esummary_params)
        summary_data = json.loads(summary_text).get("result", {})
        
        efetch_params = {
            "db": "protein",
            "id": ",".join(id_list),
            "rettype": "fasta",
            "retmode": "text"
        }
        fasta_text = http_get(EFETCH_URL, efetch_params)
        
        fasta_dict = {}
        current_id = None
        current_seq = []
        for line in fasta_text.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    fasta_dict[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            fasta_dict[current_id] = "".join(current_seq)
            
        records = []
        for uid in id_list:
            doc = summary_data.get(str(uid), {})
            accession = doc.get("accessionversion", uid)
            title = doc.get("title", "")
            seq = fasta_dict.get(accession) or fasta_dict.get(uid) or ""
            
            records.append({
                "uid": uid,
                "accession": accession,
                "title": title,
                "length": doc.get("slen", len(seq)),
                "organism": doc.get("organism", ""),
                "sequence": seq
            })

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, "a") as f:
        json.dump({"taxid": taxid, "host_filter": host_filter, "total_found": count, "records": records}, f, indent=2)
        
    print(f"[✓] Saved {len(records)} retrieved records to: {output_file}")


# =====================================================================
# STEP 2: QUALITY CONTROL (QC)
# =====================================================================
def run_qc(input_file: str, output_file: str, min_len: int = 50, max_len: int = 2000, max_ambiguous_pct: float = 1.0):
    print(f"[*] Running Quality Control on {input_file}...")
    with open(input_file, "r") as f:
        data = json.load(f)
        
    records = data.get("records", [])
    passed_records = []
    seen_sequences = set()
    
    stats = {
        "total": len(records),
        "failed_length": 0,
        "failed_ambiguity": 0,
        "failed_duplicate": 0,
        "passed": 0
    }
    
    for rec in records:
        seq = rec.get("sequence", "").upper()
        if not seq:
            continue
            
        # 1. Length check
        if len(seq) < min_len or len(seq) > max_len:
            stats["failed_length"] += 1
            continue
            
        # 2. Ambiguity check (non-standard amino acids: X, Z, B, *)
        ambiguous_count = sum(1 for aa in seq if aa not in "ACDEFGHIKLMNPQRSTVWY")
        ambiguous_pct = (ambiguous_count / len(seq)) * 100.0
        if ambiguous_pct > max_ambiguous_pct:
            stats["failed_ambiguity"] += 1
            continue
            
        # 3. Deduplication check
        if seq in seen_sequences:
            stats["failed_duplicate"] += 1
            continue
            
        seen_sequences.add(seq)
        passed_records.append(rec)
        
    stats["passed"] = len(passed_records)
    
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, "w") as f:
        json.dump({"qc_stats": stats, "records": passed_records}, f, indent=2)
        
    print(f"[✓] QC Complete. Passed {stats['passed']} / {stats['total']} records. Saved to: {output_file}")


# =====================================================================
# STEP 3: ANNOTATION & CLASSIFICATION
# =====================================================================
def annotate_proteins(input_file: str, output_dir: str):
    print(f"[*] Annotating and classifying proteins from {input_file}...")
    with open(input_file, "r") as f:
        data = json.load(f)
        
    records = data.get("records", [])
    
    # Target protein rules (Regex on title/name)
    rules = {
        "S": re.compile(r"(spike|surface glycoprotein|s protein|spike glycoprotein)", re.IGNORECASE),
        "M": re.compile(r"(\bmembrane\b|m protein|matrix protein|membrane glycoprotein)", re.IGNORECASE),
        "N": re.compile(r"(nucleocapsid|n protein|nucleoprotein)", re.IGNORECASE),
        "E": re.compile(r"(\benvelope\b|e protein|small membrane protein)", re.IGNORECASE)
    }
    
    classified = {
        "S": [],
        "M": [],
        "N": [],
        "E": [],
        "unassigned": []
    }
    
    for rec in records:
        title = rec.get("title", "")
        matched = False
        for p_class, regex in rules.items():
            if regex.search(title):
                classified[p_class].append(rec)
                matched = True
                break
        if not matched:
            classified["unassigned"].append(rec)
            
    os.makedirs(output_dir, exist_ok=True)
    
    summary = {}
    for p_class, recs in classified.items():
        summary[p_class] = len(recs)
        # Write JSON file
        json_path = os.path.join(output_dir, f"proteins_{p_class}.json")
        with open(json_path, "w") as f:
            json.dump(recs, f, indent=2)
            
        # Write FASTA file if records exist
        fasta_path = os.path.join(output_dir, f"proteins_{p_class}.fasta")
        with open(fasta_path, "w") as f:
            for r in recs:
                f.write(f">{r['accession']} {r['title']}\n{r['sequence']}\n")
                
    print(f"[✓] Annotation Complete. Summary: {summary}")
    print(f"[✓] Class-split FASTA and JSON files saved in: {output_dir}")


# =====================================================================
# STEP 4: MSA & CONSERVATIVE REGION EXTRACTION
# =====================================================================
def needleman_wunsch(seq1: str, seq2: str, match=2, mismatch=-1, gap=-2) -> Tuple[str, str]:
    """Pairwise Needleman-Wunsch global alignment."""
    n, m = len(seq1), len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1): score[i][0] = i * gap
    for j in range(m + 1): score[0][j] = j * gap
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score[i][j] = max(
                score[i-1][j-1] + s,
                score[i-1][j] + gap,
                score[i][j-1] + gap
            )
            
    align1, align2 = [], []
    i, j = n, m
    while i > 0 and j > 0:
        s = match if seq1[i-1] == seq2[j-1] else mismatch
        if score[i][j] == score[i-1][j-1] + s:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1; j -= 1
        elif score[i][j] == score[i-1][j] + gap:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    while i > 0:
        align1.append(seq1[i-1]); align2.append('-'); i -= 1
    while j > 0:
        align1.append('-'); align2.append(seq2[j-1]); j -= 1
        
    return "".join(reversed(align1)), "".join(reversed(align2))

def progressive_msa(sequences: List[str]) -> List[str]:
    """Progressive MSA aligning sequences sequentially against reference consensus."""
    if not sequences: return []
    aligned = [sequences[0]]
    for s in sequences[1:]:
        # Align new sequence against first sequence as profile anchor
        a_ref, a_new = needleman_wunsch(aligned[0], s)
        # Propagate gaps to all previously aligned sequences
        new_aligned = []
        for seq in aligned:
            seq_gapped = []
            ref_idx = 0
            for char in a_ref:
                if char == '-':
                    seq_gapped.append('-')
                else:
                    seq_gapped.append(seq[ref_idx])
                    ref_idx += 1
            new_aligned.append("".join(seq_gapped))
        new_aligned.append(a_new)
        aligned = new_aligned
    return aligned

def shannon_entropy(column: List[str]) -> Tuple[float, float]:
    """Calculates Shannon Entropy (H) and Gap ratio for an alignment column."""
    total = len(column)
    if total == 0:
        return 0.0, 1.0
        
    gaps = sum(1 for char in column if char in ("-", "."))
    gap_ratio = gaps / total
    
    residue_counts = {}
    non_gap_total = 0
    for char in column:
        if char not in ("-", "."):
            residue_counts[char] = residue_counts.get(char, 0) + 1
            non_gap_total += 1
            
    if non_gap_total == 0:
        return 0.0, gap_ratio
        
    entropy = 0.0
    for count in residue_counts.values():
        p = count / non_gap_total
        entropy -= p * math.log2(p)
        
    return entropy, gap_ratio

def run_msa_and_conservation(input_dir: str, output_dir: str, min_region_len: int = 15, max_entropy: float = 1.0):
    print(f"[*] Extracting conserved regions from classified proteins in {input_dir}...")
    os.makedirs(output_dir, exist_ok=True)
    
    classes = ["S", "M", "N", "E"]
    results = {}
    
    for p_class in classes:
        json_path = os.path.join(input_dir, f"proteins_{p_class}.json")
        if not os.path.exists(json_path):
            continue
            
        with open(json_path, "r") as f:
            recs = json.load(f)
            
        if len(recs) < 2:
            print(f"[!] Class {p_class} has fewer than 2 sequences. Skipping alignment.")
            continue
            
        sequences = [r["sequence"] for r in recs if r.get("sequence")]
        aligned_seqs = progressive_msa(sequences)
        align_len = len(aligned_seqs[0])
        
        # Calculate column-wise conservation
        column_metrics = []
        for col_idx in range(align_len):
            col = [s[col_idx] for s in aligned_seqs]
            entropy, gap_ratio = shannon_entropy(col)
            
            # Extract consensus residue
            non_gaps = [c for c in col if c not in ("-", ".")]
            consensus = max(set(non_gaps), key=non_gaps.count) if non_gaps else "-"
            
            # Require low entropy OR consensus residue present in >= 66% of sequences without excess gaps
            is_conserved = (entropy <= max_entropy) and (gap_ratio <= 0.20)
            column_metrics.append({
                "pos": col_idx + 1,
                "consensus": consensus,
                "entropy": round(entropy, 3),
                "gap_ratio": round(gap_ratio, 3),
                "is_conserved": is_conserved
            })
            
        # Find continuous conserved blocks (>= min_region_len)
        conserved_regions = []
        current_block = []
        
        for metric in column_metrics:
            if metric["is_conserved"]:
                current_block.append(metric)
            else:
                if len(current_block) >= min_region_len:
                    seq_motif = "".join(m["consensus"] for m in current_block)
                    conserved_regions.append({
                        "start": current_block[0]["pos"],
                        "end": current_block[-1]["pos"],
                        "length": len(current_block),
                        "motif": seq_motif,
                        "avg_entropy": round(sum(m["entropy"] for m in current_block) / len(current_block), 3)
                    })
                current_block = []
                
        if len(current_block) >= min_region_len:
            seq_motif = "".join(m["consensus"] for m in current_block)
            conserved_regions.append({
                "start": current_block[0]["pos"],
                "end": current_block[-1]["pos"],
                "length": len(current_block),
                "motif": seq_motif,
                "avg_entropy": round(sum(m["entropy"] for m in current_block) / len(current_block), 3)
            })
            
        results[p_class] = {
            "num_sequences": len(sequences),
            "alignment_length": align_len,
            "conserved_regions_count": len(conserved_regions),
            "conserved_regions": conserved_regions
        }
        
    out_json = os.path.join(output_dir, "conserved_regions.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
        
    print(f"[✓] Conservation extraction complete! Saved report to: {out_json}")


# =====================================================================
# CLI MAIN
# =====================================================================
def main():
    parser = argparse.ArgumentParser(description="Viral Protein Retrieval, QC, Annotation, and Conservation Pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # Fetch
    p_fetch = subparsers.add_parser("fetch", help="Fetch viral protein records from NCBI by TaxID")
    p_fetch.add_argument("--taxid", required=True, help="NCBI Taxonomy ID (e.g. 693996 for Alphacoronavirus)")
    p_fetch.add_argument("--limit", type=int, default=200, help="Max records to retrieve")
    p_fetch.add_argument("--host-filter", action="store_true", help="Filter for mammalian and pangolin hosts")
    p_fetch.add_argument("--output", required=True, help="Output JSON path")
    
    # QC
    p_qc = subparsers.add_parser("qc", help="Perform quality control filtering on retrieved proteins")
    p_qc.add_argument("--input", required=True, help="Input JSON path")
    p_qc.add_argument("--output", required=True, help="Output JSON path")
    p_qc.add_argument("--min-length", type=int, default=50, help="Minimum protein sequence length")
    p_qc.add_argument("--max-length", type=int, default=3000, help="Maximum protein sequence length")
    
    # Annotate
    p_ann = subparsers.add_parser("annotate", help="Annotate and split proteins into structural classes S, M, N, E")
    p_ann.add_argument("--input", required=True, help="Input JSON path")
    p_ann.add_argument("--output-dir", required=True, help="Output directory for annotated classes")
    
    # MSA & Conservation
    p_cons = subparsers.add_parser("msa-conservation", help="Align proteins and extract conservative regions")
    p_cons.add_argument("--input-dir", required=True, help="Directory containing annotated protein files")
    p_cons.add_argument("--output-dir", required=True, help="Output directory for conservation results")
    p_cons.add_argument("--min-region-len", type=int, default=15, help="Minimum continuous length for conserved region")
    p_cons.add_argument("--max-entropy", type=float, default=0.3, help="Max Shannon entropy for conservation threshold")
    
    args = parser.parse_args()
    
    if args.command == "fetch":
        # Retrieve the sequences from NCBI using the provided Taxonomy ID array in the CSV file
        df = pd.read_csv(args.taxid)
        for taxid in df['Taxonomy ID'].tolist():
            fetch_viral_proteins(taxid, args.limit, args.output, host_filter=args.host_filter)
    elif args.command == "qc":
        run_qc(args.input, args.output, min_len=args.min_length, max_len=args.max_length)
    elif args.command == "annotate":
        annotate_proteins(args.input, args.output_dir)
    elif args.command == "msa-conservation":
        run_msa_and_conservation(args.input_dir, args.output_dir, min_region_len=args.min_region_len, max_entropy=args.max_entropy)

if __name__ == "__main__":
    main()
