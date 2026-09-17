# conserved_regions.py

Find conserved peptide-length windows in a protein multiple sequence alignment.

Takes an existing MSA, decides which columns are conserved by plain per-column
amino-acid identity, then reports every window of a given length range that falls
entirely inside a run of conserved columns. The windows are intended as **peptide
candidates**, which is why overlapping sub-windows are emitted rather than one row
per region.

The tool does not align anything. Feed it an alignment from MAFFT, MUSCLE or
similar.

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

Turn it on as a sanity check: if a block survives weighting, it is real
conservation; if it disappears, it was an artifact of over-represented near-duplicate
entries. Databases like RefSeq are full of these for well-sampled viruses, and
exact-duplicate removal does not catch them.

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

At genus level, drop `--min-identity` to around 0.70 or group the input more narrowly.
`--columns-out` tells you which threshold is binding. `PIPELINE_INTEGRATION.md` §10.1
has the full numbers.

## Tests

```sh
python3 test_conserved_regions.py
```

33 assertions with hand-computed expected values, covering tiling arithmetic, both
gap-voting modes, the occupancy guard, ambiguity handling, degenerate columns,
threshold inclusivity, strict contiguity, and reference-coordinate mapping. Exits
non-zero on failure.
