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

### Per-sequence matrix (`--matrix-out`)

The window table's `motif` is a **consensus** and need not equal any one sequence.
This file says what each input sequence actually carries. One row per sequence, three
columns per conserved stretch, and **two header lines** — read it with
`pandas.read_csv(path, sep="\t", header=[0, 1], index_col=0)`.

```
region   	5-18          	5-18    	5-18    	21-30     	21-30   	21-30
field    	aligned       	residues	identity	aligned   	residues	identity
consensus	ACDEFGHIKLMNPQ	NA      	NA      	RSTVWYACDE	NA      	NA
sp_A     	ACDEFGHIKLMNPQ	4-17    	1.0000  	RSTVWYACDE	18-27   	1.0000
sp_C     	ACDEFSHIKLMNPQ	3-16    	0.9286  	RSTVWYACDE	19-28   	1.0000
sp_E     	ACDEFGHIKLMNPQ	4-17    	1.0000  	RST---ACDE	18-24   	1.0000
sp_F     	ACDEFGHIKLMNPQ	5-18    	1.0000  	----------	NA      	NA
```

| field | notes |
|---|---|
| header line 1 | the stretch's extent in alignment columns, repeated three times |
| header line 2 | which of the three values the column holds |
| `consensus` row | what `identity` is measured against; `residues`/`identity` are `NA` |
| `aligned` | that sequence's slice **verbatim from the alignment**, gaps included, always the width of the stretch |
| `residues` | the same span in that sequence's **own ungapped numbering**, `first-last`, or `NA` if it is all gaps there |
| `identity` | fraction of this sequence's residues in the stretch that match the consensus, 4 dp, or `NA` |

Stretch labels and ordering match the `stretch`/`stretch_start`/`stretch_end` fields
of the window table, so the files join.

`identity` uses the same denominator rule as `--gap-votes`: by default a sequence with
an internal deletion is scored over the residues it does have (`sp_E` above reads
`1.0000` — no substitutions), and with `--gap-votes` the missing positions count
against it (`0.7000`). An ambiguity code never matches the consensus.

### Per-sequence FASTA (`--fasta-out`)

The same subsequences as sequence records, **gaps stripped**, grouped by stretch:

```
>sp_C_S1 region=5-18 residues=3-16 identity=0.9286
ACDEFSHIKLMNPQ
>sp_E_S2 region=21-30 residues=18-24 identity=1.0000
RSTACDE
```

The first token is a unique ID (`<input id>_S<n>`); the rest are `key=value`. A
sequence that is all gaps in a stretch gets **no record** — the count of skipped
sequences goes to stderr, so a short file is not a silent failure.

**Records within a stretch are not all the same length.** Stripping gaps is right for
peptides you intend to order, but it means `sp_E_S2` above is 7 aa from a 10-column
stretch. If you need column correspondence — a sequence logo, say — use the `aligned`
column of the matrix instead, not this file.

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
This follows the specification and is not a bug, but on real viral alignments it does
fragment usable candidates — see §10.1.

The obvious fix is a `--bridge N` option merging runs separated by at most N failing
columns. **Do not implement that as stated.** On the betacoronavirus spike alignment
in §10.1, a simulated `--bridge 2` at `--min-identity 0.95` turns 0 windows into 196 —
but the merged regions contain columns at 25–45% identity:

```
region(cols)   len   worst column identity inside
 1238-1263      26   0.450
 1286-1295      10   0.450
 1358-1384      27   0.250
 1389-1402      14   0.450
```

A window advertised as 95%-conserved that contains a 25% column is not a peptide
candidate, and nothing downstream can tell the difference. If bridging is added,
bridged columns need their own identity floor (a second, lower threshold) and the
window row needs to report it.

**`motif` is the consensus, not a real strain.** At ≥95% identity these are nearly
always identical but can differ at up to 5% of positions. If you are ordering
peptides and need the exact residues of a particular isolate, take them from that
sequence using the `--reference` coordinates.

### 10.1 Behaviour on real data

Checked 2026-09-17 against RefSeq structural proteins (S, E, M, N) for two
coronavirus genera, length/ambiguity/duplicate filtered, aligned with MAFFT `--auto`.

| group | seqs | columns | windows at defaults |
|---|---|---|---|
| *Alphacoronavirus* S / E / M / N | 31 / 29 / 28 / 28 | 1736 / 89 / 263 / 552 | 0 / 0 / 0 / 0 |
| *Betacoronavirus* S / E / M / N | 20 / 14 / 19 / 18 | 1627 / 90 / 231 / 526 | 0 / 0 / 0 / 0 |

**Zero windows at the defaults, in all eight groups.** This is correct behaviour, not
a failure. Conserved columns do exist — 158 in *Alphacoronavirus* S, 182 in
*Betacoronavirus* S — but they never form a run of 10. The longest contiguous run at
95% identity is 5 columns. Genus-wide 95% identity across 10+ consecutive positions
essentially does not occur.

Longest contiguous conserved run / windows emitted, by identity threshold:

```
group           95%       90%       80%       70%       60%
alphacov S      5/0       5/0       8/0      12/6      14/30
alphacov E      1/0       1/0       1/0       2/0       3/0
alphacov M      3/0       3/0       6/0       8/0      16/42
alphacov N      3/0       5/0       6/0       6/0      12/6
betacov  S      5/0       6/0      10/1      11/3      12/15
betacov  E      1/0       1/0       1/0       3/0       3/0
betacov  M      3/0       3/0       8/0       8/0      13/10
betacov  N      3/0       7/0       7/0       7/0      12/6
```

**The defaults are scoped to within-species or within-subgenus comparisons.** If you
are working at genus level, either relax `--min-identity` to 0.70 or below, or group
the input more narrowly. Use `--columns-out` to see which of the two thresholds is
binding before changing either.

**It does recover real biology.** The one region found in *Betacoronavirus* S at
`--min-identity 0.70` is alignment columns 1324–1334, consensus `AQIDRLINGR`. Against
`--reference YP_009724390.1` that is SARS-CoV-2 spike residues 991–1000 — in context
`ILSRLDKVEAE·VQIDRLITGR·LQSLQTYVTQQLIRA`, the HR1 / central helix of S2 and a
well-characterised broadly-neutralising target. Out of 1627 columns spanning a whole
genus, that is the region it selected.

**E is not productive at this scale.** 4–5 conserved columns out of ~90, never more
than 3 consecutive even at 60% identity. Short and divergent; do not expect candidates
from it above species level.

### 10.2 Species level, and why you must use `--weighting henikoff`

Same checks against MERS-CoV (taxid 1335626), GenBank rather than RefSeq, 250 entries
fetched per protein, filtered on length and ambiguity, **exact duplicates removed**:

| | unique seqs | columns | conserved cols | windows at defaults |
|---|---|---|---|---|
| S | 82 | 1470 | 1056 | 1143 |
| E | 9 | 82 | 58 | 1 |
| M | 26 | 339 | 75 | 1 |
| N | 49 | 443 | 340 | 446 |

The defaults behave as intended here — the same flags that return nothing across a
genus return 1143 spike windows within a species. That is the scoping conclusion from
§10.1 confirmed from the other direction.

**The important result is what weighting does to those numbers:**

```
        weighting=none   weighting=henikoff
S       1143 windows  ->  159 windows    (-86%)
N        446 windows  ->  109 windows    (-76%)
M          1 window   ->    0 windows    (gone entirely)
E          1 window   ->    1 window     (unchanged)
```

All exact duplicates had already been removed — 168 of the 250 fetched spike entries.
**That was not close to sufficient.** The 82 remaining *unique* sequences still carry
enough near-duplicate structure to inflate the window count roughly sevenfold. M is
the clearest case: its single conserved 10-mer disappears completely under weighting,
so it was never conservation, just one strain sampled repeatedly with a few
substitutions. E, with 9 sequences and no redundancy left, is identical under both
settings and acts as the control.

Sampling bias in GenBank is not random — outbreak strains are sequenced hundreds of
times — so it manufactures apparent conservation in exactly the regions least likely
to be broadly useful. **For any GenBank-derived set, run `--weighting henikoff`, and
treat the unweighted count as an upper bound rather than a result.**

The 159 surviving spike windows collapse to 11 regions, against
`--reference sp|K9N5Q8.1|SPIKE_MERS1` (SwissProt canonical, strain HCoV-EMC/2012):

| region | MERS-S | motif | annotation |
|---|---|---|---|
| S1 | 401–410 | `RLVFTNCNYN` | RBD core (367–588) |
| S2 | 799–808 | `IQKVTVDCKQ` | S2, upstream of the S2′ site |
| S3 | 899–908 | `TIADPGYMQG` | fusion peptide (888–910) |
| S4 | 929–938 | `VAGYKVLPPL` | post-fusion-peptide |
| S5 | 984–993 | `GITQQVLSEN` | HR1 start (984–1104) |
| S6 | 997–1006 | `IANKFNQALG` | HR1 |
| S7 | 1029–1038 | `NAQALSKLAS` | HR1 |
| S8 | 1066–1075 | `QIDRLINGRL` | central helix |
| S9 | 1078–1087 | `LNAFVAQQLV` | central helix |
| S10 | 1276–1285 | `LNESYIDLKE` | HR2 (1246–1295) |
| S11 | 1327–1336 | `CMGKLKCNRC` | cysteine-rich cytoplasmic tail |

Eight of eleven fall in S2, the conserved fusion machinery; S1, which carries the
variable antigenic surface, contributes one. That distribution is what a working
conservation scan should produce.

Note S8 against §10.1: the genus-level *Betacoronavirus* scan independently returned
`AQIDRLINGR` at SARS-CoV-2 S 991–1000. These are the same central-helix site,
recovered from different input at a different taxonomic scale, sharing `QIDRLI`
exactly. Two independent paths to the same residues is a useful end-to-end check on
scoring and on the reference-coordinate mapping.

### 10.3 What per-sequence output exposes that the aggregate tables hide

Re-running the Henikoff-weighted MERS-CoV spike set (82 sequences, 11 stretches, 889
subsequence records) with `--matrix-out` surfaced three things invisible in the
window table. None is a tool bug; all three are things to check before ordering
anything.

**1. 15 of 889 records do not match their own stretch's consensus.** The window table
reports `RLVFTNCNYN` for the RBD stretch at columns 439–448. `YEV46216.1` actually
carries `RLIFTNCNYN` there (`identity=0.9000`) — verified against its raw ungapped
sequence at residues 401–410. A peptide synthesised from the consensus does not cover
that strain. Most stretches are clean, but `S7` (columns 1097–1116) has five
imperfect records and `S1`, `S2`, `S3`, `S6` have one or two each.

**2. Four PDB entries in the set are construct fragments, not spike proteins.** In
`S11` (columns 1414–1425, cysteine-rich cytoplasmic tail) the consensus is
`CMGKLKCNRCCD`, but:

```
pdb|9DKK|A   GLNDI          identity=0.0000
pdb|7YMX|C   DNSAD          identity=0.2000
pdb|22FY|B   LEVL           identity=0.2500
pdb|5W9I|J   LEVL           identity=0.2500
```

A record at 0.0000 identity inside a stretch that passed a 95% identity filter is the
signature of a structure-derived construct — expression tags and linkers aligned into
a region the construct does not actually span. Filter `pdb|` accessions at retrieval,
or accept that they contribute junk to the ragged ends.

**3. Henikoff weight is heavily concentrated, and that is worth looking at.** Four of
82 sequences carry **39% of the total vote**; `USL83011.1` alone carries 17.3%. We
checked whether this was an artifact of low-quality records — it is not: those four
are full-length and carry no ambiguity codes, so the weighting is behaving exactly as
designed, rewarding genuine divergence. But it does mean the "95%" in a weighted run
is a consensus of a handful of divergent sequences plus a large near-identical block,
and §10.2's recommendation to always weight GenBank sets should be read alongside
this. The `identity` column is how you see it; the window table cannot show it.

### 10.4 Note for the retrieval step

Two NCBI E-utilities behaviours cost real time here, both upstream of this program but
worth recording. `[Protein Name]` is not a valid search field — queries using it
return zero hits silently rather than erroring. And a multi-word `[All Fields]` term
is tokenised, so `membrane protein[All Fields]` matches `membrane AND protein` and
pulled a 76 aa envelope protein into an M-protein set. Length bounds caught it, which
is an argument for keeping per-protein length QC even when the query looks specific.

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
python3 test_conserved_regions.py     # 55 assertions, exits non-zero on failure
```

Covers tiling arithmetic, both gap-voting modes, the occupancy guard, ambiguity
handling, degenerate columns, threshold inclusivity, strict contiguity,
reference-coordinate mapping, per-sequence identity (including internal deletions
under both gap-voting modes) and the layout of both per-sequence output files. Worth running in CI — the assertions have hand-computed
expected values, so a failure means real behaviour changed.
