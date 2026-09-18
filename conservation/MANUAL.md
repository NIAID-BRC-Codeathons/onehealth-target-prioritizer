# conserved_regions.py — reference manual

Complete description of every setting and every output.

For a short introduction see `README.md`. For the guarantees other pipeline
components can rely on, and for what the tool does on real coronavirus data, see
`PIPELINE_INTEGRATION.md`.

---

## Contents

1. [Synopsis](#1-synopsis)
2. [Input](#2-input)
3. [How conservation is decided](#3-how-conservation-is-decided)
4. [Settings](#4-settings)
5. [Outputs](#5-outputs)
6. [Exit codes and diagnostics](#6-exit-codes-and-diagnostics)
7. [A worked example, end to end](#7-a-worked-example-end-to-end)
8. [Recipes](#8-recipes)
9. [Performance](#9-performance)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Synopsis

```
conserved_regions.py ALIGNMENT [options]
```

Reads a protein multiple sequence alignment. Decides which alignment columns are
conserved, by plain per-column amino-acid identity. Finds maximal runs of
consecutive conserved columns ("stretches"). Emits every sub-window of a given
length range inside each stretch, as peptide candidates.

It does **not** align, retrieve, group, filter, rank or annotate, and it never
touches the network. Pure Python, NumPy and Biopython; deterministic, so identical
input gives byte-identical output on any machine.

**Requirements:** Python 3.11+, NumPy, Biopython.

---

## 2. Input

### File and format

The single positional argument is the alignment file. Format is guessed from the
extension, or forced with `--format`.

| extension | format passed to Biopython |
|---|---|
| `.fa` `.fasta` `.fas` `.faa` `.afa` `.mfa` | `fasta` |
| `.aln` `.clustal` `.clw` | `clustal` |
| `.sto` `.stk` `.stockholm` | `stockholm` |
| `.phy` `.phylip` | `phylip-relaxed` |
| `.nex` `.nexus` | `nexus` |
| anything else, including no extension | `fasta` |

The fallback is `fasta`, so an unrecognised extension is not an error — it is an
assumption. If the file is not FASTA, say so with `--format`.

`--format` accepts any format name Biopython's `AlignIO` supports, not only the ones
listed above.

### What the input must satisfy

| requirement | if violated |
|---|---|
| at least 2 sequences | exit 1, `alignment needs at least 2 sequences, got N` |
| all sequences the same length | exit 1, `Sequences must all be the same length` |
| parseable as the chosen format | exit 1, the parser's own message |

**The input must already be aligned.** This tool does not align anything, and it
cannot tell a bad alignment from a good one — it will happily report conservation
from a misaligned file. Use MAFFT, MUSCLE or an equivalent.

Sequence identifiers are taken from Biopython's `record.id`, which for FASTA is the
first whitespace-delimited token of the header line.

### Alphabet handling

Every character is mapped to one of 22 internal codes.

| input | treated as | notes |
|---|---|---|
| `ACDEFGHIKLMNPQRSTVWY` | that amino acid | the 20 standard residues |
| lowercase `acdefghiklmnpqrstvwy` | the same amino acid | case is not meaningful; soft-masking is ignored |
| `-` `.` `~` | gap | all three are equivalent |
| `X` `B` `Z` `U` `O` `*` | ambiguity | see below |
| any other symbol | ambiguity | there is no "invalid character" error |

**Ambiguity codes occupy a position but can never win it.** They count toward
occupancy, and they count in the identity denominator, but they cannot become the
consensus residue. A column of 6 `X` and 4 `A` scores identity 0.4, not 1.0, and its
consensus is `A`.

The alternative — ignoring ambiguous symbols entirely — would let a column with two
known residues and 98 unknowns score 100% identity. That is the wrong answer for a
tool whose output gets synthesised.

---

## 3. How conservation is decided

A column is **conserved** when both tests pass:

```
identity   >= --min-identity     the most common amino acid's share of the column
occupancy  >= --min-occupancy    the share of sequences with any residue there
```

Both comparisons are **inclusive**: a column at exactly 0.95 passes `--min-identity
0.95`.

### The identity denominator, and `--gap-votes`

| mode | denominator | column of 90 `A` + 10 gaps |
|---|---|---|
| `--no-gap-votes` (default) | non-gap symbols only | 90/90 = **1.0000**, conserved at 95% |
| `--gap-votes` | all sequences | 90/100 = **0.9000**, not conserved at 95% |

A gap is never counted as "the same amino acid" in either mode. It can only ever
affect the denominator.

Under `--gap-votes`, `--min-occupancy` becomes largely redundant: a column cannot
reach 95% identity unless it is at least 95% occupied anyway. The occupancy guard
matters in the default mode.

### Why the occupancy guard exists

Without it, a column where only 3 of 200 sequences have a residue scores 100%
identity and is called perfectly conserved. Ragged alignment ends produce a lot of
those. The default 0.5 floor discards them. `--min-occupancy 0` turns the guard off.

### Sequence weighting

With `--weighting henikoff`, every sequence gets a weight from Henikoff & Henikoff
(1994) position-based weighting, and all counts above become weighted sums instead
of plain counts. Weights are normalised to sum to the number of sequences, so a
uniform weighting reproduces the unweighted result exactly.

Columns with only one symbol type contribute nothing, since an invariant column says
nothing about redundancy.

The effect is to down-weight over-represented near-identical sequences. This matters
enormously on database-derived sets — see §4.6 and `PIPELINE_INTEGRATION.md` §10.2.

### Ties

When two amino acids are equally frequent, the one with the lower index in
`ACDEFGHIKLMNPQRSTVWY` wins. This keeps output reproducible across runs and machines.

It is worth stating because the obvious idiom, `max(set(col), key=col.count)`, is
*not* reproducible: set iteration order for strings varies with `PYTHONHASHSEED`.

### From columns to windows

1. Mark every column conserved or not.
2. Find **maximal runs** of consecutive conserved columns. These are *stretches*,
   labelled `S1`, `S2`, … in column order.
3. Discard stretches shorter than `--min-length`.
4. **Tile** each surviving stretch: emit every sub-window of length `--min-length`
   through `--max-length` that fits inside it.

Runs are **strictly contiguous**. One failing column splits a stretch in two, and
each fragment is then length-filtered independently. This is deliberate; see §10 of
`PIPELINE_INTEGRATION.md` for why bridging was evaluated and rejected.

---

## 4. Settings

### Summary

| option | type | default | effect |
|---|---|---|---|
| `ALIGNMENT` | path | *required* | the input MSA |
| `--format` | string | guessed from extension | input format |
| `--min-identity` | fraction | `0.95` | identity threshold for a column |
| `--min-occupancy` | fraction | `0.5` | occupancy threshold for a column |
| `--gap-votes` / `--no-gap-votes` | flag | `--no-gap-votes` | gaps in the identity denominator |
| `--weighting` | `none` \| `henikoff` | `none` | sequence weighting |
| `--min-length` | int | `10` | shortest window, and minimum stretch length |
| `--max-length` | int | `15` | longest window |
| `--reference` | sequence id | none | add coordinates in one sequence's numbering |
| `--windows-out` | path | stdout | where the window table goes |
| `--columns-out` | path | not written | per-column profile |
| `--matrix-out` | path | not written | per-sequence × per-stretch table |
| `--fasta-out` | path | not written | per-sequence subsequences as FASTA |
| `-h`, `--help` | flag | | usage and exit 0 |

### 4.1 `--format FORMAT`

Overrides extension-based guessing. Use it for files with an unusual or absent
extension, since the fallback is `fasta` and a wrong guess produces a parse error
rather than a helpful one.

### 4.2 `--min-identity F`

The most common amino acid must hold at least this share of the column.

**Accepts either `0.95` or `95`.** Any value above 1.0 is divided by 100. The result
must land in `[0, 1]`; `--min-identity 150` is a usage error (exit 2).

Lowering this is the single biggest lever on output volume, and the effect is
strongly non-linear — see §9.

The default assumes input from **one species or subgenus**. Across whole genera it is
too strict and will legitimately return nothing; see §10.

### 4.3 `--min-occupancy F`

A column is only eligible if at least this share of sequences has a non-gap symbol
there. Same `0.5` / `50` parsing as `--min-identity`.

Set `0` to disable the guard entirely. Raise it toward `1.0` to demand that a region
be present in nearly every sequence, which is usually what you want if the windows
are going to be synthesised.

### 4.4 `--gap-votes` / `--no-gap-votes`

Whether gaps count in the identity denominator. Default is `--no-gap-votes`.

This flag also governs the `identity` figures in the `--matrix-out` and `--fasta-out`
files, so the whole output set stays internally consistent — see §5.3.

Use `--gap-votes` when a region must be *present*, not merely *unvaried where present*.

### 4.5 `--weighting {none,henikoff}`

`none` (default) means every sequence counts once, so "95%" means 95% of the
sequences in the file. That is the number you want in a methods section.

`henikoff` applies position-based weighting, which down-weights redundant
near-identical sequences.

### 4.6 A note on when to use weighting

**Turn weighting on for anything pulled from GenBank.** Removing exact duplicates is
not sufficient. On 82 *unique* MERS-CoV spike sequences — every exact duplicate
already stripped — weighting cut the output from 1143 windows to 159, and the single
window found for the M protein disappeared entirely, because it was one strain
sequenced repeatedly rather than conservation.

Treat an unweighted count from a database-derived set as an upper bound.

The counterweight: weighting concentrates influence as hard as redundancy dilutes it.
On that same set, 4 of 82 sequences carry 39% of the total weight and one carries
17.3%. That is the algorithm working as designed — rewarding genuine divergence — but
it means the consensus is driven by a handful of sequences. Use `--matrix-out` to see
which. Full numbers in `PIPELINE_INTEGRATION.md` §10.2 and §10.3.

### 4.7 `--min-length N` / `--max-length N`

`--min-length` does double duty: it is both the shortest window emitted **and** the
minimum stretch length, since a stretch shorter than the minimum window cannot
contain one.

`--max-length` must be greater than or equal to `--min-length`, or the run is a usage
error (exit 2). `--min-length` must be at least 1.

For a single window length, set both: `--min-length 12 --max-length 12`.

### 4.8 `--reference SEQ_ID`

Adds two columns to the window table, and one to the column profile, giving positions
in that sequence's **own ungapped residue numbering** rather than alignment columns.

The identifier may be a **unique prefix** of a sequence id. An unknown id and an
ambiguous prefix are both exit 1, and the ambiguity message lists up to five matches.

This affects reporting only. It never changes which columns are conserved.

Note that `--matrix-out` and `--fasta-out` report own-numbering coordinates for
**every** sequence, so `--reference` is not needed for those.

### 4.9 Output destinations

`--windows-out` redirects the window table from stdout to a file. The other three are
purely additive: nothing is written unless you ask.

All four can be used together in a single run, which is cheaper than four runs since
scoring happens once.

---

## 5. Outputs

All tables are tab-separated with a header. **All coordinates are 1-based and
inclusive**, so a span `5-18` covers 14 columns, and `length` is `end - start + 1`.

### 5.1 Window table — stdout, or `--windows-out`

One row per window. This is the primary output.

| field | type | notes |
|---|---|---|
| `id` | `W<n>` | sequential across the whole run, from `W1` |
| `stretch` | `S<n>` | which maximal run this window came from |
| `stretch_start`, `stretch_end` | int | that run's full extent, alignment columns |
| `start`, `end` | int | this window's extent, alignment columns |
| `length` | int | `end - start + 1` |
| `<ref>_start`, `<ref>_end` | int or `NA` | **only with `--reference`**; the column name contains the id you passed |
| `min_identity` | float, 4 dp | lowest column identity in the window |
| `mean_identity` | float, 4 dp | mean column identity across the window |
| `min_occupancy` | float, 4 dp | lowest column occupancy in the window |
| `motif` | string | consensus residues, length `length` |

Two things to watch:

**The reference columns are conditionally present and dynamically named.** If you
parse by position this will break. Parse by header name.

**`motif` is a consensus, not any real strain's residues.** At ≥95% identity these are
nearly always the same, but they can differ at up to 5% of positions. If you need the
exact residues of a particular sequence, use `--matrix-out` or `--fasta-out`.

### 5.2 Column profile — `--columns-out`

One row per alignment column, whether conserved or not. The diagnostic view: use it
when a region you expected did not appear, to see which threshold was binding.

| field | type | notes |
|---|---|---|
| `column` | int | 1-based alignment column |
| `ref_pos` | int or `NA` | **only with `--reference`**; `NA` where the reference is gapped |
| `consensus` | char | most common amino acid, or `-` if the column has none |
| `identity` | float, 4 dp | as defined in §3, reflecting `--gap-votes` and `--weighting` |
| `occupancy` | float, 4 dp | non-gap share |
| `conserved` | `1` / `0` | both thresholds met |

A column of nothing but gaps and/or ambiguity codes has no consensus and is reported
as `-`, rather than an arbitrary amino acid.

### 5.3 Per-sequence matrix — `--matrix-out`

The two tables above are aggregates. This one says what each individual sequence
actually carries, which is what you need in order to pick a real strain to order a
peptide from.

**Layout:** one row per input sequence, three columns per stretch, and **two header
lines**. Total width is `1 + 3 × (number of stretches)`.

```
region     5-18            5-18      5-18      21-30       21-30     21-30
field      aligned         residues  identity  aligned     residues  identity
consensus  ACDEFGHIKLMNPQ  NA        NA        RSTVWYACDE  NA        NA
sp_A       ACDEFGHIKLMNPQ  4-17      1.0000    RSTVWYACDE  18-27     1.0000
sp_B       ACDEFGHIKLMNPQ  5-18      1.0000    RSTVWYACDE  19-28     1.0000
sp_C       ACDEFSHIKLMNPQ  3-16      0.9286    RSTVWYACDE  19-28     1.0000
sp_D       ACDEFGHIKLMNPQ  5-18      1.0000    RSTVWYASDE  19-28     0.9000
sp_E       ACDEFGHIKLMNPQ  4-17      1.0000    RST---ACDE  18-24     1.0000
sp_F       ACDEFGHIKLMNPQ  5-18      1.0000    ----------  NA        NA
```

(Shown space-aligned for readability; the file itself is tab-separated.)

| row / field | meaning |
|---|---|
| header line 1 | the stretch's extent in alignment columns, repeated three times |
| header line 2 | which of the three values the column holds |
| `consensus` row | what `identity` is measured against; its `residues` and `identity` cells are `NA` |
| `aligned` | that sequence's slice **verbatim from the alignment**, gaps included, always exactly the width of the stretch |
| `residues` | the same span in that sequence's **own ungapped numbering**, as `first-last`, or `NA` if it is all gaps there |
| `identity` | fraction of this sequence's residues in the stretch matching the consensus, 4 dp, or `NA` |

Read it with:

```python
pandas.read_csv(path, sep="\t", header=[0, 1], index_col=0)
```

Stretch labels and ordering match the window table's `stretch`, `stretch_start` and
`stretch_end`, so the files join.

**`identity` follows `--gap-votes`.** By default a sequence with an internal deletion
is scored over the residues it does have — `sp_E` above reads `1.0000`, meaning no
substitutions. Under `--gap-votes` the missing positions count against it and it reads
`0.7000`. One says "no substitutions", the other says "don't order this one". Pick
whichever your downstream step needs; the `aligned` cell shows the gaps either way, so
nothing is hidden.

An ambiguity code never matches the consensus, so it scores as a mismatch.

If there are no stretches, the file contains the header lines and the sequence names,
with no data columns.

### 5.4 Per-sequence FASTA — `--fasta-out`

The same subsequences as sequence records, **gaps stripped**, grouped by stretch.

```
>sp_C_S1 region=5-18 residues=3-16 identity=0.9286
ACDEFSHIKLMNPQ
>sp_E_S2 region=21-30 residues=18-24 identity=1.0000
RSTACDE
```

Description line grammar:

```
>{sequence id}_S{stretch number} region={start}-{end} residues={first}-{last} identity={fraction}
```

The first token is a unique record id; the rest are `key=value` pairs. `region` is in
alignment columns, `residues` in that sequence's own numbering.

Records are ordered by stretch, then by input order within each stretch.

**A sequence that is all gaps in a stretch gets no record.** The number skipped is
reported on stderr, so a short file is not a silent failure.

**Records within a stretch are not all the same length.** Stripping gaps is right for
peptides you intend to order — but it means `sp_E_S2` above is 7 aa from a 10-column
stretch. If you need column correspondence, for a sequence logo or a sub-alignment,
use the `aligned` column of `--matrix-out` instead.

If there are no stretches, the file is empty (0 bytes).

### 5.5 stderr summary

Two lines, always, plus a third when `--fasta-out` is used:

```
[*] 82 sequences x 1470 columns | identity>=95% occupancy>=50% gaps=ignored weighting=henikoff
[*] 663 conserved columns -> 11 stretches >= 10 aa -> 159 windows of 10-15 aa
[*] 889 subsequence records across 11 stretches (13 skipped, all gaps in that stretch)
```

The first line echoes the settings in force, which makes a log self-documenting. The
second gives the funnel: conserved columns, then stretches surviving the length
filter, then windows. Reading it left to right tells you which step lost your region.

stderr is safe to discard. **stdout carries the window table and nothing else**, so
`> windows.tsv` cleanly separates them.

---

## 6. Exit codes and diagnostics

| code | meaning |
|---|---|
| `0` | ran successfully — **including "found zero windows"** |
| `1` | input error: missing file, unparseable, ragged, unknown or ambiguous `--reference` |
| `2` | usage error from argparse: unknown flag, `--max-length < --min-length`, out-of-range fraction |

**Zero windows is not an error.** It is a result, and often the correct one. The exit
code stays 0 and the window table is written with its header and no rows. If you need
to branch on "found something", count rows — do not test the exit status.

---

## 7. A worked example, end to end

Six sequences, thirty columns. Small enough to verify every number by hand.

```
                        1111111111222222222233
               1234567890123456789012345678901
    sp_A       MKT-ACDEFGHIKLMNPQ--RSTVWYACDE
    sp_B       MKTLACDEFGHIKLMNPQ--RSTVWYACDE
    sp_C       M--QACDEFSHIKLMNPQGGRSTVWYACDE     S at column 10
    sp_D       MKTWACDEFGHIKLMNPQ--RSTVWYASDE     S at column 28
    sp_E       MKT-ACDEFGHIKLMNPQ--RST---ACDE     internal deletion
    sp_F       MKTRACDEFGHIKLMNPQ------------     truncated
```

Run at 80% identity, because with six sequences a single mismatch is 83% and the
default 95% would reject every interesting case:

```sh
./conserved_regions.py mock.fasta --min-identity 0.80 --reference sp_C \
    --columns-out cols.tsv --matrix-out matrix.tsv --fasta-out sub.faa
```

### stderr

```
[*] 6 sequences x 30 columns | identity>=80% occupancy>=50% gaps=ignored weighting=none
[*] 27 conserved columns -> 2 stretches >= 10 aa -> 16 windows of 10-15 aa
[*] 11 subsequence records across 2 stretches (1 skipped, all gaps in that stretch)
```

### Column profile (first eight rows)

```
column  ref_pos  consensus  identity  occupancy  conserved
1       1        M          1.0000    1.0000     1
2       NA       K          1.0000    0.8333     1
3       NA       T          1.0000    0.8333     1
4       2        L          0.2500    0.6667     0
5       3        A          1.0000    1.0000     1
6       4        C          1.0000    1.0000     1
7       5        D          1.0000    1.0000     1
8       6        E          1.0000    1.0000     1
```

Column 2 shows the occupancy denominator at work: five sequences have `K`, one has a
gap, so identity is 5/5 but occupancy is 5/6.

Column 4 shows tie-breaking. The residues are `L`, `Q`, `W`, `R`, one each. All four
tie at 0.25, and `L` wins because it has the lowest index in
`ACDEFGHIKLMNPQRSTVWY`.

Columns 1–3 form a conserved run of length 3. It never becomes a stretch, because 3
is below `--min-length 10`.

### Window table (16 rows)

```
id   stretch  stretch_start  stretch_end  start  end  length  sp_C_start  sp_C_end  min_identity  mean_identity  min_occupancy  motif
W1   S1       5              18           5      14   10      3           12        0.8333        0.9833         1.0000         ACDEFGHIKL
W2   S1       5              18           5      15   11      3           13        0.8333        0.9848         1.0000         ACDEFGHIKLM
W3   S1       5              18           5      16   12      3           14        0.8333        0.9861         1.0000         ACDEFGHIKLMN
...
W16  S2       21             30           21     30   10      19          28        0.8000        0.9800         0.6667         RSTVWYACDE
```

`min_identity 0.8333` in `S1` is column 10, where `sp_C` has `S` against five `G`.
The reference columns are named `sp_C_start` and `sp_C_end` because that is the id
passed to `--reference`.

### Per-sequence matrix

```
region     5-18            5-18      5-18      21-30       21-30     21-30
field      aligned         residues  identity  aligned     residues  identity
consensus  ACDEFGHIKLMNPQ  NA        NA        RSTVWYACDE  NA        NA
sp_A       ACDEFGHIKLMNPQ  4-17      1.0000    RSTVWYACDE  18-27     1.0000
sp_B       ACDEFGHIKLMNPQ  5-18      1.0000    RSTVWYACDE  19-28     1.0000
sp_C       ACDEFSHIKLMNPQ  3-16      0.9286    RSTVWYACDE  19-28     1.0000
sp_D       ACDEFGHIKLMNPQ  5-18      1.0000    RSTVWYASDE  19-28     0.9000
sp_E       ACDEFGHIKLMNPQ  4-17      1.0000    RST---ACDE  18-24     1.0000
sp_F       ACDEFGHIKLMNPQ  5-18      1.0000    ----------  NA        NA
```

Every residue span is checkable against the alignment above. `sp_C` begins `M--Q`, so
its numbering runs three lower than `sp_B`'s, and the table shows exactly that.

`sp_C` scores 13/14 = 0.9286 in `S1`. `sp_D` scores 9/10 = 0.9000 in `S2`. `sp_E`
scores 7/7 = 1.0000 in `S2` under the default denominator, and would score 0.7000
under `--gap-votes`.

### Per-sequence FASTA

```
>sp_A_S1 region=5-18 residues=4-17 identity=1.0000
ACDEFGHIKLMNPQ
...
>sp_D_S2 region=21-30 residues=19-28 identity=0.9000
RSTVWYASDE
>sp_E_S2 region=21-30 residues=18-24 identity=1.0000
RSTACDE
```

Eleven records, not twelve: `sp_F` is all gaps in `S2`, so it is skipped and counted
on stderr. `sp_E_S2` is **7 aa from a 10-column stretch**, because its gaps are gone.

---

## 8. Recipes

**Collapse the tiled windows back to distinct regions**

```sh
./conserved_regions.py aln.fasta | awk 'NR>1 {print $2"\t"$3"-"$4}' | sort -u
```

**One window length only**

```sh
./conserved_regions.py aln.fasta --min-length 12 --max-length 12
```

**Everything, in one pass**

```sh
./conserved_regions.py aln.fasta --weighting henikoff \
    --windows-out windows.tsv --columns-out columns.tsv \
    --matrix-out matrix.tsv --fasta-out subseqs.faa
```

**Find sequences that do not carry their own region's consensus**

```sh
grep '^>' subseqs.faa | grep -v 'identity=1.0000'
```

Each hit is a strain the consensus peptide would not cover.

**Demand a region be present, not merely unvaried**

```sh
./conserved_regions.py aln.fasta --gap-votes --min-occupancy 0.9
```

**Find out which threshold is blocking a region you expected**

```sh
./conserved_regions.py aln.fasta --columns-out cols.tsv
awk -F'\t' 'NR>1 && $1>=400 && $1<=450' cols.tsv
```

Compare the `identity` and `occupancy` columns against your thresholds.

---

## 9. Performance

Scoring is vectorised over columns. The MERS-CoV spike set (82 sequences ×
1470 columns) completes in about 0.15 s including all four output files, on a laptop.

Two things dominate at scale:

**Henikoff weighting iterates over columns in Python**, so its cost grows linearly
with alignment width. It is still cheap at typical protein sizes.

**Window count grows quadratically with stretch length**, and this is what actually
bites. Every sub-window of every stretch is emitted:

| stretch length | windows at 10–15 aa |
|---:|---:|
| 9 | 0 |
| 10 | 1 |
| 12 | 6 |
| 15 | 21 |
| 20 | 51 |
| 40 | 171 |
| 100 | 531 |

A well-conserved 1500-column alignment can produce well over 1000 rows. That is
intended — the rows are orderable peptide candidates, and selection happens
downstream — but plan for the volume, and note that `--matrix-out` is one column-pair
per *stretch*, not per window, so it stays small.

---

## 10. Troubleshooting

**Zero windows, and I expected some.**

Most likely correct, and the most common cause is that the input spans too broad a
taxonomic range. Across all eight coronavirus genus-level groups tested, the defaults
return zero, because the longest run of 95%-identical columns anywhere is 5. The same
flags return 1143 spike windows within a single species.

Check `--columns-out` before assuming breakage. If conserved columns exist but no
stretch reaches `--min-length`, that is a real biological result, not a bug. Group the
input more narrowly, or lower `--min-identity` — but see the next entry.

**Lowering the threshold produced lots of windows. Is that a result?**

Not by itself. Check what is inside them. Lowering a threshold until output appears is
not the same as finding conservation, and a window is only as good as its worst
column. `min_identity` in the window table is the number to look at.

**A "conserved" region vanishes when I add `--weighting henikoff`.**

Then it was an artifact of redundant near-identical sequences, and the weighted answer
is the right one. Exact-duplicate removal does not catch this. See §4.6.

**A sequence scores 0.0000 identity inside a conserved stretch.**

Look at that record. On real data this turned out to be structure-derived construct
entries — expression tags and linkers aligned into a region the construct does not
actually span. Consider filtering `pdb|` accessions at retrieval.

**My FASTA has fewer records than sequences × stretches.**

Expected: sequences that are all gaps in a stretch are skipped, and the count is on
stderr. Ragged alignment ends and partial records are the usual cause.

**The reference columns moved, or my parser broke.**

They are conditionally present and dynamically named after the `--reference` id.
Parse by header name, not column position.

**Everything is reported as conserved.**

Check for `--min-occupancy 0` combined with the default `--no-gap-votes`. That
combination will call a column of one residue and 199 gaps perfectly conserved.

---

## Tests

```sh
python3 test_conserved_regions.py
```

55 assertions with hand-computed expected values, covering tiling arithmetic, both
gap-voting modes, the occupancy guard, ambiguity handling, degenerate columns,
threshold inclusivity, strict contiguity, reference-coordinate mapping, per-sequence
identity including internal deletions under both gap-voting modes, and the layout of
both per-sequence output files. Exits non-zero on failure.
