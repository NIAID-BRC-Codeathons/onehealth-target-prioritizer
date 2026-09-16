# conserved_regions.py — what it is and how to build around it

A note for anyone wiring this into a larger pipeline. `README.md` covers day-to-day
usage; this document covers the **contract**: what the program guarantees, what it
expects, and what it deliberately leaves to you.

---

## 1. What it does, in one paragraph

Given a protein multiple sequence alignment, it labels each alignment column
conserved or not by plain amino-acid identity, finds maximal runs of consecutive
conserved columns, and emits every fixed-length window that fits inside those runs.
The windows are **peptide candidates** — that is the reason overlapping sub-windows
are emitted instead of one row per region. Output is a TSV.

## 2. Where it sits

```
  sequence retrieval  ->  QC / dedup  ->  grouping  ->  ALIGNER  ->  conserved_regions.py  ->  downstream
  (NCBI, local FASTA)     (length,        (by gene,     (MAFFT,      (this program)            (peptide
                           ambiguity,      protein,      MUSCLE,                                 selection,
                           duplicates)     clade...)     ...)                                    ranking,
                                                                                                 ordering)
```

It occupies exactly one box. Everything to its left and right is yours to build.

## 3. What it does NOT do — read this first

These are scope boundaries, not missing features. Assuming otherwise is the most
likely way to misuse it.

| It does not | Consequence for you |
|---|---|
| **Align anything** | You must supply an already-aligned MSA. Use MAFFT or MUSCLE. |
| **Group sequences** | It scores whatever you hand it as one alignment. Splitting by gene/protein/clade happens upstream. |
| **Touch the network** | No NCBI calls, no downloads. Fully offline and reproducible. |
| **Filter or QC sequences** | Garbage sequences in the alignment become garbage columns. QC upstream. |
| **Rank or deduplicate windows** | Overlapping windows are emitted in full. Selection and ranking are downstream. |
| **Write temp files or need external binaries** | Pure Python + NumPy + Biopython. |
| **Handle nucleotides** | Protein only. DNA/RNA would be read as ambiguity codes and score nonsense. |

## 4. Input contract

- **Formats:** FASTA, Clustal, Stockholm, PHYLIP (relaxed), NEXUS. Guessed from the
  file extension, overridable with `--format`.
- **Must be aligned.** All sequences the same length, or it exits with an error. It
  will not pad or trim for you.
- **At least 2 sequences.**
- **Gap characters:** `-`, `.`, `~` are all treated as gaps.
- **Case-insensitive.** Lowercase residues are treated as residues, not as masked
  positions. If your upstream soft-masks with lowercase, that information is lost
  here — hard-mask to `X` instead if it matters.
- **Non-standard symbols** (`X B Z U O *`, and anything non-ASCII) are treated as
  ambiguity codes: they occupy a position but can never be the consensus residue.
- **Sequence IDs** are Biopython's `record.id`, i.e. the FASTA header up to the first
  whitespace. That is what `--reference` matches against (unique prefixes are
  accepted).

## 5. Output contract

### Process behaviour

| | |
|---|---|
| **stdout** | the window TSV, and nothing else — safe to pipe |
| **stderr** | a two-line human summary — safe to discard |
| **exit 0** | ran successfully (including "found zero windows") |
| **exit 1** | input error: missing file, unparseable, ragged, unknown `--reference` |
| **exit 2** | usage error from argparse (bad flags, `--max-length < --min-length`) |

Zero windows is **not** an error and does not change the exit code. Check the row
count, not the exit status, if you need to branch on "found something".

### Window table

Tab-separated, one header line, one row per window.

| field | type | notes |
|---|---|---|
| `id` | `W<n>` | sequential over the whole run, starting at `W1` |
| `stretch` | `S<n>` | which maximal conserved run this window came from |
| `stretch_start`, `stretch_end` | int | that run's full extent |
| `start`, `end` | int | this window's extent |
| `length` | int | `end - start + 1` |
| `<ref>_start`, `<ref>_end` | int or `NA` | **only present with `--reference`**; the column name contains the sequence ID you passed |
| `min_identity`, `mean_identity` | float, 4 dp | over the window's columns |
| `min_occupancy` | float, 4 dp | lowest occupancy in the window |
| `motif` | string | consensus residues, length == `length` |

**Coordinates are 1-based and inclusive, in alignment columns.** Reference
coordinates are 1-based ungapped residue positions in that sequence. Both ends are
inclusive, so `length` is `end - start + 1`, not `end - start`.

Note that the reference columns are **conditionally present and dynamically named**.
If you parse by position this will break; parse by header name.

### Column profile (`--columns-out`)

One row per alignment column: `column`, `ref_pos` (only with `--reference`),
`consensus`, `identity`, `occupancy`, `conserved` (`1`/`0`). This is the diagnostic
view — use it when a region you expected did not appear, to see whether the cause was
low identity or low occupancy.

## 6. The conservation rule

Precise statement, so you can describe it in a methods section.

For each column, let **W** be the total sequence weight (with default weighting, the
number of sequences):

- **top** = the largest weight held by any one of the 20 standard amino acids
- **occupied** = total weight of all non-gap symbols, ambiguity codes included
- **occupancy** = occupied / W
- **identity** = top / occupied  (default), or top / W  (with `--gap-votes`)

A column is **conserved** when `identity >= --min-identity` **and**
`occupancy >= --min-occupancy`. Both comparisons are inclusive: exactly 95% passes at
a 95% threshold.

Gaps never count as "the same amino acid" in either mode — they can only affect the
denominator. Ambiguity codes count toward `occupied` but can never supply `top`.

A window's reported `min_identity` is the minimum over its columns, so every column
in a window is guaranteed to meet the threshold.

## 7. Parameters

| flag | default | effect |
|---|---|---|
| `--min-identity` | `0.95` | column identity threshold; accepts `95` or `0.95` |
| `--min-occupancy` | `0.5` | minimum non-gap fraction for a column to be eligible |
| `--gap-votes` / `--no-gap-votes` | gaps ignored | whether gaps enter the identity denominator |
| `--weighting` | `none` | `henikoff` down-weights redundant near-identical sequences |
| `--min-length` / `--max-length` | `10` / `15` | window length range |
| `--reference` | none | adds coordinates in one sequence's own numbering |
| `--windows-out` / `--columns-out` | stdout / none | output destinations |
| `--format` | guessed | input format override |

### Predicting output size

For a conserved run of length **L**, the number of windows is

```
sum over w from min_length to min(max_length, L) of (L - w + 1)
```

With the 10–15 defaults: a run of 10 gives 1 window, 12 gives 6, 15 gives 21, 22
gives 63, 40 gives 171, 100 gives 531. A long conserved region produces a lot of
rows. Set `--min-length` and `--max-length` equal for one window length only.

To collapse back to distinct biological regions:

```sh
conserved_regions.py aln.fasta | awk 'NR>1 {print $2"\t"$3"-"$4}' | sort -u
```

## 8. Two ways to call it

### As a subprocess

```python
import csv, io, subprocess

proc = subprocess.run(
    ["python3", "conserved_regions.py", "aln.fasta",
     "--reference", "NP_828851.1", "--min-identity", "95"],
    capture_output=True, text=True, check=True,
)
windows = list(csv.DictReader(io.StringIO(proc.stdout), delimiter="\t"))
for w in windows:
    print(w["motif"], w["start"], w["end"], w["NP_828851.1_start"])
```

`check=True` raises on exit 1 and 2. `proc.stderr` holds the summary.

### As a module

Import it if you want the arrays rather than a TSV. This is the full pipeline in
eight lines, and it is exactly what the CLI does:

```python
import conserved_regions as cr

aln     = cr.load_alignment("aln.fasta")                 # Alignment(ids, codes)
weights = cr.WEIGHTINGS["none"](aln.codes)               # or "henikoff"
scores  = cr.score_columns(aln.codes, weights,
                           min_identity=0.95, min_occupancy=0.5, gap_votes=False)
runs    = cr.maximal_runs(scores.conserved, min_length=10)

for start, end in runs:                                  # 0-based, inclusive
    for ws, we in cr.tile(start, end, 10, 15):
        print(ws + 1, we + 1, scores.motif(ws, we))      # +1 to match the TSV
```

**Module-level indices are 0-based; the TSV is 1-based.** The CLI adds 1 on output.

Useful objects:

- `Alignment.codes` — `(n_seq, n_col)` int8. `0–19` are the standard amino acids in
  `STANDARD_AA` order, `20` is any ambiguity code, `21` is a gap.
- `Alignment.index_of(seq_id)` — row index, accepts a unique prefix, raises `KeyError`.
- `ColumnScores.identity`, `.occupancy`, `.conserved` — NumPy arrays of length `n_col`.
- `ColumnScores.consensus` — consensus string of length `n_col`.
- `reference_map(codes_row)` — column index to 1-based ungapped position, `0` at gaps.

## 9. Performance

200 sequences × 1300 columns runs in ~0.13 s. Memory is `n_seq × n_col` bytes for the
alignment plus a small `n_col × 22` float matrix. A 5000 × 3000 alignment is roughly
15 MB. Neither time nor memory should constrain pipeline design at realistic sizes.

Runs are **deterministic**: identical input and flags give byte-identical output
across runs and machines. Ties between equally frequent amino acids are broken by
position in `ACDEFGHIKLMNPQRSTVWY`. There is no randomness and no hash-order
dependence anywhere.

## 10. Known limitations

**No bridging.** Runs must be strictly contiguous, so one variable column splits a
stretch in two and each fragment is length-filtered independently. A single
hypervariable position inside an otherwise conserved 30-mer costs the whole region.
This follows the specification and is not a bug, but on real viral alignments it will
fragment usable candidates. The fix, if needed, is a `--bridge N` option merging runs
separated by at most N failing columns.

**`motif` is the consensus, not a real strain.** At ≥95% identity these are nearly
always identical but can differ at up to 5% of positions. If you are ordering
peptides and need the exact residues of a particular isolate, take them from that
sequence using the `--reference` coordinates.

**Untested on real data.** Every check so far has used synthetic alignments. The
defaults are reasonable but unvalidated against real viral alignments, where
redundancy and ragged ends behave differently. Use `--columns-out` to check the
thresholds are landing sensibly before trusting the first real run.

## 11. Upstream requirement: use a real aligner

The alignment quality sets a ceiling on everything here — this program cannot detect
conservation that the aligner destroyed.

Do not hand-roll a progressive aligner for this pipeline. A star-alignment approach
that aligns each new sequence against a single reference and propagates gaps
back through the previously aligned set is easy to get subtly wrong: gap columns
introduced by one sequence must be distinguished from gap columns already present in
the reference, and if they are not, indices drift and residues are silently dropped
from sequences aligned earlier in the process. That kind of bug produces no
exception — the output looks like a normal alignment, just a wrong one.

For the alignment step, call an established aligner (MAFFT or MUSCLE) rather than
implementing one. Retrieval, QC and annotation are separate, tractable problems and
can reasonably be custom code; alignment is where the general-purpose tools already
handle the hard cases (guide trees, affine gap penalties, iterative refinement) that
a first-pass implementation is likely to get wrong.

## 12. Tests

```sh
python3 test_conserved_regions.py     # 33 assertions, exits non-zero on failure
```

Covers tiling arithmetic, both gap-voting modes, the occupancy guard, ambiguity
handling, degenerate columns, threshold inclusivity, strict contiguity, and
reference-coordinate mapping. Worth running in CI — the assertions have hand-computed
expected values, so a failure means real behaviour changed.
