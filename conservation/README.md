# conserved_regions.py

Find conserved peptide-length windows in a protein multiple sequence alignment.

Takes an existing MSA, decides which columns are conserved by plain per-column
amino-acid identity, then reports every window of a given length range that falls
entirely inside a run of conserved columns. The windows are intended as **peptide
candidates**, which is why overlapping sub-windows are emitted rather than one row
per region.

The tool does not align anything. Feed it an alignment from MAFFT, MUSCLE or
similar.

Full reference for every setting and every output field: **[`MANUAL.md`](MANUAL.md)**.
Guarantees for neighbouring pipeline components, and results on real data:
[`PIPELINE_INTEGRATION.md`](PIPELINE_INTEGRATION.md).

## Requirements

Python 3.11+, NumPy, Biopython.

```sh
pip install numpy biopython
```

## Quick start

```sh
# defaults: identity >= 95%, occupancy >= 50%, gaps ignored, windows 10-15 aa
./conserved_regions.py aln.fasta

# report coordinates in one strain's own numbering, and keep the column profile
./conserved_regions.py aln.fasta --reference NP_828851.1 --columns-out profile.tsv

# stricter: make gaps count against a column, and relax identity to 90%
./conserved_regions.py aln.fasta --gap-votes --min-identity 90

# what each individual sequence carries in each conserved stretch
./conserved_regions.py aln.fasta --matrix-out regions.tsv --fasta-out regions.faa
```

Input format is guessed from the extension (`.fasta` `.fa` `.faa` `.aln` `.sto`
`.phy` `.nex`), or set it with `--format`. The window table goes to stdout; a one
line summary goes to stderr, so `> windows.tsv` keeps them separate.

## How conservation is defined

A column is **conserved** when both hold:

| test | meaning | default |
|---|---|---|
| identity | the most common amino acid's share of the column | ≥ 95% |
| occupancy | fraction of sequences with any residue (non-gap) there | ≥ 50% |

`--min-identity` and `--min-occupancy` accept either form: `95` and `0.95` mean the
same thing.

### Gap voting

`--gap-votes` decides what goes in the identity denominator.

| mode | denominator | column of 90 A + 10 gaps |
|---|---|---|
| `--no-gap-votes` (default) | non-gap symbols only | 90/90 = **100%**, conserved |
| `--gap-votes` | all sequences | 90/100 = **90%**, not conserved at 95% |

A gap is never counted as "the same amino acid" in either mode — it can only ever
affect the denominator.

Note that `--gap-votes` makes `--min-occupancy` largely redundant: with gaps voting,
a column cannot reach 95% identity unless it is at least 95% occupied anyway. The
guard matters in the default mode.

### Occupancy guard

Without it, a column where only 3 of 200 sequences have a residue would score 100%
identity and be called perfectly conserved. Ragged alignment ends produce a lot of
those. The default 50% floor discards them; `--min-occupancy 0` turns the guard off.

### Sequence weighting

`--weighting henikoff` applies Henikoff & Henikoff (1994) position-based weights,
which down-weight redundant near-identical sequences. It is **off by default** so
that "95%" means 95% of the sequences in the file, which is what you want in a
methods section.

Off by default, but **turn it on for anything pulled from GenBank.** If a block
survives weighting it is real conservation; if it disappears, it was an artifact of
over-represented near-duplicate entries, and exact-duplicate removal does not catch
those.

This is not a marginal correction. On 82 *unique* MERS-CoV spike sequences — every
exact duplicate already stripped — weighting cuts the output from 1143 windows to 159,
and the single window found for the M protein disappears entirely. Outbreak strains
get sequenced hundreds of times, and that sampling bias manufactures apparent
conservation. Treat the unweighted count as an upper bound. Full numbers in
`PIPELINE_INTEGRATION.md` §10.2.

## Windows

Maximal runs of consecutive conserved columns are found first, then every window of
length `--min-length`..`--max-length` (default 10–15) inside each run is emitted.

This is deliberately verbose — the point is to enumerate peptide candidates:

| conserved run | windows emitted |
|---:|---:|
| 9 columns | 0 (shorter than the minimum) |
| 10 | 1 |
| 12 | 6 |
| 15 | 21 |
| 20 | 51 |
| 22 | 63 |
| 40 | 171 |
| 100 | 531 |

Every row carries `stretch`, `stretch_start` and `stretch_end`, so collapsing back to
distinct biological regions is a `sort -u` away:

```sh
./conserved_regions.py aln.fasta | awk 'NR>1 {print $2"\t"$3"-"$4}' | sort -u
```

For a single window length, set both bounds: `--min-length 12 --max-length 12`.

## Output

### Window table (stdout, or `--windows-out`)

| column | meaning |
|---|---|
| `id` | window id, `W1`, `W2`, … |
| `stretch` | which maximal conserved run this window came from, `S1`, `S2`, … |
| `stretch_start`, `stretch_end` | that run's full extent, 1-based alignment columns |
| `start`, `end`, `length` | this window, 1-based alignment columns |
| `<ref>_start`, `<ref>_end` | same span in the reference's own residue numbering (only with `--reference`; `NA` if the reference is all gaps there) |
| `min_identity`, `mean_identity` | identity across the window's columns |
| `min_occupancy` | lowest occupancy in the window |
| `motif` | consensus residues across the window |

### Column profile (`--columns-out`)

One row per alignment column: `column`, `ref_pos` (with `--reference`), `consensus`,
`identity`, `occupancy`, `conserved`. Useful for plotting the profile or for
debugging why a block did not survive.

### Per-sequence matrix (`--matrix-out`)

The two tables above are aggregates. `motif` is a consensus and need not match any one
sequence — so if you are picking a strain to order a peptide from, it does not tell you
what that strain actually has. This file does.

One row per input sequence, three columns per conserved stretch, two header lines:

```
region   	5-18          	5-18    	5-18    	21-30     	21-30   	21-30
field    	aligned       	residues	identity	aligned   	residues	identity
consensus	ACDEFGHIKLMNPQ	NA      	NA      	RSTVWYACDE	NA      	NA
sp_A     	ACDEFGHIKLMNPQ	4-17    	1.0000  	RSTVWYACDE	18-27   	1.0000
sp_C     	ACDEFSHIKLMNPQ	3-16    	0.9286  	RSTVWYACDE	19-28   	1.0000
sp_E     	ACDEFGHIKLMNPQ	4-17    	1.0000  	RST---ACDE	18-24   	1.0000
sp_F     	ACDEFGHIKLMNPQ	5-18    	1.0000  	----------	NA      	NA
```

| | |
|---|---|
| `aligned` | the slice verbatim from the alignment, gaps included |
| `residues` | the same span in that sequence's own ungapped numbering, or `NA` |
| `identity` | fraction of this sequence's residues that match the consensus, or `NA` |

Load it with `pandas.read_csv(path, sep="\t", header=[0, 1], index_col=0)`.

`identity` follows the `--gap-votes` rule like everything else: by default `sp_E` scores
`1.0000` because its seven residues all match, and under `--gap-votes` it scores
`0.7000` because three positions are missing. The first says "no substitutions", the
second says "don't order this one" — pick the one your downstream step needs.

### Per-sequence FASTA (`--fasta-out`)

The same subsequences, **gaps stripped**, grouped by stretch:

```
>sp_C_S1 region=5-18 residues=3-16 identity=0.9286
ACDEFSHIKLMNPQ
>sp_E_S2 region=21-30 residues=18-24 identity=1.0000
RSTACDE
```

A sequence that is all gaps in a stretch gets no record; the skipped count goes to
stderr. Because gaps are stripped, records within a stretch can differ in length —
`sp_E_S2` is 7 aa from a 10-column stretch. That is what you want for ordering, and
wrong for a sequence logo; use the matrix's `aligned` column for anything that needs
column correspondence.

## Design decisions

These are deliberate. Each is easy to change if it does not suit your data.

**Ambiguity codes occupy a position but never win it.** `X`, `B`, `Z`, `U`, `O` and
`*` count toward occupancy and toward the identity denominator, but cannot be the
consensus residue. So a column of 6 X + 4 A scores 0.4, not 1.0, and its consensus is
`A`. The alternative — skipping ambiguous symbols entirely — would let a column with
two known residues and 98 unknowns score 100%.

**A column of nothing but gaps and/or ambiguity codes has no consensus** and is
reported as `-` rather than an arbitrary amino acid.

**`motif` is the consensus sequence, not any real strain's residues.** At ≥95%
identity these are nearly always identical, but they can differ at up to 5% of
positions. If you are ordering peptides and need the exact residues of a particular
strain, take them from that sequence using the `<ref>_start`/`<ref>_end` coordinates.

**Ties are broken deterministically.** When two amino acids are equally frequent the
lower index in `ACDEFGHIKLMNPQRSTVWY` wins, so output is reproducible across runs and
machines. (This is worth stating because the obvious idiom,
`max(set(col), key=col.count)`, is *not* reproducible — set iteration order for
strings varies with `PYTHONHASHSEED`.)

## Known limitation: no bridging

Runs must be strictly contiguous. A single variable column splits a stretch in two,
and each fragment is then length-filtered independently. On real viral alignments
this will fragment otherwise-usable peptide candidates — one hypervariable position
inside an otherwise conserved 30-mer costs you the whole region.

This follows the specification ("continuous stretches of conserved columns") and is
not a bug. The obvious fix is a `--bridge N` option merging runs separated by at most
N failing columns, but on real data that turns out to be a bad trade: bridging a
betacoronavirus spike alignment with `N=2` at 95% identity produces windows containing
columns at 25–45% identity. A "95%-conserved" peptide with a 25% position in the
middle is not a usable candidate. See `PIPELINE_INTEGRATION.md` §10 before building
it.

## Choosing thresholds

The 95% / 10–15 aa defaults assume sequences from **one species or subgenus**. They
are too strict for a whole genus: across all eight coronavirus genus-level groups
tested (S, E, M, N for *Alphacoronavirus* and *Betacoronavirus*), the defaults return
zero windows, because the longest run of 95%-identical columns anywhere is 5.

Within a species they work as intended: the same flags that return nothing across
*Betacoronavirus* return 1143 spike windows on MERS-CoV.

At genus level, drop `--min-identity` to around 0.70 or group the input more narrowly.
`--columns-out` tells you which threshold is binding. `PIPELINE_INTEGRATION.md` §10.1
and §10.2 have the full numbers for both scales.

## Tests

```sh
python3 test_conserved_regions.py
```

55 assertions with hand-computed expected values, covering tiling arithmetic, both
gap-voting modes, the occupancy guard, ambiguity handling, degenerate columns,
threshold inclusivity, strict contiguity, reference-coordinate mapping, per-sequence
identity (including internal deletions under both gap-voting modes) and the layout of
both per-sequence output files. Exits non-zero on failure.
