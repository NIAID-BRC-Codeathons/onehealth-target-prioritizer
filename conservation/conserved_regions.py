#!/usr/bin/env python3
"""Find conserved peptide-length windows in a protein multiple sequence alignment.

Rule
----
A column is CONSERVED when the most common amino acid occupies at least
--min-identity of the column, and the column is occupied (non-gap) in at least
--min-occupancy of the sequences.

    gaps not voting (default)  identity = top_aa / (non-gap symbols)
    gaps voting (--gap-votes)  identity = top_aa / (all sequences)

A gap is never counted as "the same amino acid", and neither is an ambiguity
code (X/B/Z/U/O/*): ambiguous symbols occupy a position, so they count in the
denominator, but they can never win the column.

Maximal runs of consecutive conserved columns are then tiled into every window
of length --min-length .. --max-length that fits inside them, which is what you
want when the windows are peptide candidates rather than domain annotations.

Usage
-----
    ./conserved_regions.py aln.fasta
    ./conserved_regions.py aln.fasta --min-identity 90 --gap-votes
    ./conserved_regions.py aln.fasta --reference NP_828851.1 \
        --weighting henikoff --columns-out profile.tsv
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from functools import cached_property
from typing import Iterator, TextIO

import numpy as np
from Bio import AlignIO

# --------------------------------------------------------------------------
# Alphabet
# --------------------------------------------------------------------------
# Columns are coded as small integers. 0-19 are the standard amino acids, 20 is
# any other residue symbol (X, B, Z, U, O, *), 21 is a gap.

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
OTHER_CODE = 20
GAP_CODE = 21
N_CODES = 22
GAP_CHARS = "-.~"

_CODE_OF = np.full(128, OTHER_CODE, dtype=np.int8)
for _i, _aa in enumerate(STANDARD_AA):
    _CODE_OF[ord(_aa)] = _i
    _CODE_OF[ord(_aa.lower())] = _i
for _g in GAP_CHARS:
    _CODE_OF[ord(_g)] = GAP_CODE

_FORMAT_BY_EXT = {
    "fa": "fasta", "fasta": "fasta", "fas": "fasta", "faa": "fasta",
    "afa": "fasta", "mfa": "fasta",
    "aln": "clustal", "clustal": "clustal", "clw": "clustal",
    "sto": "stockholm", "stk": "stockholm", "stockholm": "stockholm",
    "phy": "phylip-relaxed", "phylip": "phylip-relaxed",
    "nex": "nexus", "nexus": "nexus",
}


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------
@dataclass
class Alignment:
    ids: list[str]
    codes: np.ndarray  # (n_seq, n_col) int8

    @property
    def n_seq(self) -> int:
        return self.codes.shape[0]

    @property
    def n_col(self) -> int:
        return self.codes.shape[1]

    def index_of(self, seq_id: str) -> int:
        """Resolve a sequence id, allowing a unique prefix match."""
        if seq_id in self.ids:
            return self.ids.index(seq_id)
        hits = [i for i, name in enumerate(self.ids) if name.startswith(seq_id)]
        if len(hits) == 1:
            return hits[0]
        if not hits:
            raise KeyError(f"sequence {seq_id!r} not found in alignment")
        raise KeyError(f"sequence {seq_id!r} is ambiguous, matches "
                       + ", ".join(self.ids[i] for i in hits[:5]))


def load_alignment(path: str, fmt: str | None = None) -> Alignment:
    fmt = fmt or _FORMAT_BY_EXT.get(path.rsplit(".", 1)[-1].lower(), "fasta")
    records = list(AlignIO.read(path, fmt))
    if len(records) < 2:
        raise ValueError(f"alignment needs at least 2 sequences, got {len(records)}")

    widths = {len(r.seq) for r in records}
    if len(widths) != 1:
        raise ValueError(f"sequences have differing lengths {sorted(widths)} "
                         "-- input must be aligned")

    raw = np.array([bytearray(str(r.seq), "ascii") for r in records], dtype=np.uint8)
    raw[raw > 127] = ord("X")
    return Alignment(ids=[r.id for r in records], codes=_CODE_OF[raw])


# --------------------------------------------------------------------------
# Sequence weighting (optional; off by default so "95%" means 95% of sequences)
# --------------------------------------------------------------------------
def henikoff_weights(codes: np.ndarray) -> np.ndarray:
    """Henikoff & Henikoff (1994) position-based weights, normalised to sum to n_seq.

    Down-weights over-represented near-identical sequences, which otherwise
    dominate every column and inflate apparent conservation.
    """
    n_seq, n_col = codes.shape
    w = np.zeros(n_seq, dtype=np.float64)
    for j in range(n_col):
        col = codes[:, j]
        counts = np.bincount(col, minlength=N_CODES)
        n_types = int((counts > 0).sum())
        if n_types > 1:  # invariant columns say nothing about redundancy
            w += 1.0 / (n_types * counts[col])
    if not w.any():
        return np.ones(n_seq, dtype=np.float64)
    return w * (n_seq / w.sum())


WEIGHTINGS = {
    "none": lambda codes: np.ones(codes.shape[0], dtype=np.float64),
    "henikoff": henikoff_weights,
}


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------
@dataclass
class ColumnScores:
    """Per-column identity, occupancy and consensus.

    identity   fraction of the column held by its most common amino acid
    occupancy  fraction of sequences with any residue (non-gap) in the column
    top_code   index into STANDARD_AA of that most common amino acid
    has_aa     whether any standard amino acid is present at all
    conserved  identity and occupancy both met their thresholds
    """

    identity: np.ndarray
    occupancy: np.ndarray
    top_code: np.ndarray
    has_aa: np.ndarray
    conserved: np.ndarray

    @cached_property
    def consensus(self) -> str:
        # A column of nothing but gaps and/or ambiguity codes has no consensus.
        return "".join(STANDARD_AA[c] if ok else "-"
                       for c, ok in zip(self.top_code, self.has_aa))

    def motif(self, start: int, end: int) -> str:
        return self.consensus[start : end + 1]


def score_columns(codes: np.ndarray, weights: np.ndarray, *, min_identity: float,
                  min_occupancy: float, gap_votes: bool) -> ColumnScores:
    n_col = codes.shape[1]
    counts = np.zeros((n_col, N_CODES), dtype=np.float64)
    for code in range(N_CODES):
        mask = codes == code
        if mask.any():
            counts[:, code] = weights @ mask

    total = float(weights.sum())
    # Only the 20 standard amino acids can win a column; ambiguity codes and
    # gaps cannot, but ambiguity codes still occupy the position.
    aa_counts = counts[:, :OTHER_CODE]
    top_count = aa_counts.max(axis=1)
    top_code = aa_counts.argmax(axis=1)  # first-max, so ties are deterministic

    occupied = counts[:, :GAP_CODE].sum(axis=1)
    occupancy = occupied / total

    denom = np.full(n_col, total) if gap_votes else occupied
    identity = np.divide(top_count, denom, out=np.zeros(n_col), where=denom > 0)

    return ColumnScores(
        identity=identity,
        occupancy=occupancy,
        top_code=top_code,
        has_aa=top_count > 0,
        conserved=(identity >= min_identity) & (occupancy >= min_occupancy),
    )


# --------------------------------------------------------------------------
# Region calling
# --------------------------------------------------------------------------
def maximal_runs(conserved: np.ndarray, min_length: int) -> list[tuple[int, int]]:
    """Maximal runs of consecutive True, as inclusive 0-based (start, end) pairs."""
    idx = np.flatnonzero(conserved)
    if idx.size == 0:
        return []
    runs: list[list[int]] = [[int(idx[0]), int(idx[0])]]
    for i in idx[1:]:
        i = int(i)
        if i == runs[-1][1] + 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return [(s, e) for s, e in runs if e - s + 1 >= min_length]


def tile(start: int, end: int, min_length: int, max_length: int) -> Iterator[tuple[int, int]]:
    """Every sub-window of length min_length..max_length inside [start, end].

    Yielded in start-major order, shortest first at each start.
    """
    for s in range(start, end - min_length + 2):
        longest = min(max_length, end - s + 1)
        for w in range(min_length, longest + 1):
            yield s, s + w - 1


def reference_map(codes_row: np.ndarray) -> np.ndarray:
    """Column index -> 1-based ungapped position, 0 where the reference is gapped."""
    non_gap = codes_row != GAP_CODE
    return np.where(non_gap, np.cumsum(non_gap), 0)


def reference_span(ref_pos: np.ndarray, start: int, end: int) -> tuple[str, str]:
    present = ref_pos[start : end + 1][ref_pos[start : end + 1] > 0]
    if present.size == 0:
        return "NA", "NA"
    return str(int(present[0])), str(int(present[-1]))


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def write_windows(out: TextIO, runs, scores: ColumnScores, min_length: int,
                  max_length: int, ref_pos: np.ndarray | None,
                  ref_name: str | None) -> int:
    header = ["id", "stretch", "stretch_start", "stretch_end", "start", "end", "length"]
    if ref_pos is not None:
        header += [f"{ref_name}_start", f"{ref_name}_end"]
    header += ["min_identity", "mean_identity", "min_occupancy", "motif"]
    print("\t".join(header), file=out)

    n = 0
    for s_num, (r_start, r_end) in enumerate(runs, 1):
        for w_start, w_end in tile(r_start, r_end, min_length, max_length):
            n += 1
            ident = scores.identity[w_start : w_end + 1]
            occ = scores.occupancy[w_start : w_end + 1]
            row = [f"W{n}", f"S{s_num}", str(r_start + 1), str(r_end + 1),
                   str(w_start + 1), str(w_end + 1), str(w_end - w_start + 1)]
            if ref_pos is not None:
                row += list(reference_span(ref_pos, w_start, w_end))
            row += [f"{ident.min():.4f}", f"{ident.mean():.4f}", f"{occ.min():.4f}",
                    scores.motif(w_start, w_end)]
            print("\t".join(row), file=out)
    return n


def write_columns(out: TextIO, scores: ColumnScores, ref_pos: np.ndarray | None) -> None:
    header = ["column", "consensus", "identity", "occupancy", "conserved"]
    if ref_pos is not None:
        header.insert(1, "ref_pos")
    print("\t".join(header), file=out)

    consensus = scores.consensus
    for j in range(scores.identity.size):
        row = [str(j + 1)]
        if ref_pos is not None:
            row.append(str(int(ref_pos[j])) if ref_pos[j] else "NA")
        row += [consensus[j], f"{scores.identity[j]:.4f}",
                f"{scores.occupancy[j]:.4f}", "1" if scores.conserved[j] else "0"]
        print("\t".join(row), file=out)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def fraction(text: str) -> float:
    """Accept either 0.95 or 95 for a percentage-style option."""
    value = float(text)
    if value > 1.0:
        value /= 100.0
    if not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError(f"{text!r} is not a fraction or percentage")
    return value


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Find conserved peptide-length windows in a protein MSA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("alignment", help="input MSA file")
    p.add_argument("--format", help="input format (default: guessed from extension)")

    g = p.add_argument_group("column conservation")
    g.add_argument("--min-identity", type=fraction, default=0.95, metavar="F",
                   help="fraction of the column that must share one amino acid "
                        "(accepts 0.95 or 95)")
    g.add_argument("--min-occupancy", type=fraction, default=0.5, metavar="F",
                   help="fraction of sequences that must be non-gap for a column "
                        "to be eligible at all")
    g.add_argument("--gap-votes", action=argparse.BooleanOptionalAction, default=False,
                   help="count gaps in the identity denominator")
    g.add_argument("--weighting", default="none", choices=sorted(WEIGHTINGS),
                   help="sequence weighting; 'henikoff' offsets redundant sequences "
                        "but makes the percentage a weighted one")

    g = p.add_argument_group("windows")
    g.add_argument("--min-length", type=int, default=10, help="shortest window")
    g.add_argument("--max-length", type=int, default=15, help="longest window")

    g = p.add_argument_group("output")
    g.add_argument("--reference", metavar="SEQ_ID",
                   help="also report coordinates in this sequence's own numbering")
    g.add_argument("--windows-out", metavar="PATH",
                   help="write the window table here (default: stdout)")
    g.add_argument("--columns-out", metavar="PATH",
                   help="also write the full per-column profile here")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.min_length < 1:
        parser.error("--min-length must be >= 1")
    if args.max_length < args.min_length:
        parser.error("--max-length must be >= --min-length")

    try:
        aln = load_alignment(args.alignment, args.format)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    ref_pos = ref_name = None
    if args.reference:
        try:
            ref_pos = reference_map(aln.codes[aln.index_of(args.reference)])
        except KeyError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        ref_name = args.reference

    weights = WEIGHTINGS[args.weighting](aln.codes)
    scores = score_columns(
        aln.codes, weights,
        min_identity=args.min_identity,
        min_occupancy=args.min_occupancy,
        gap_votes=args.gap_votes,
    )
    runs = maximal_runs(scores.conserved, args.min_length)

    out = open(args.windows_out, "w") if args.windows_out else sys.stdout
    try:
        n_windows = write_windows(out, runs, scores, args.min_length,
                                  args.max_length, ref_pos, ref_name)
    finally:
        if args.windows_out:
            out.close()

    print(
        f"[*] {aln.n_seq} sequences x {aln.n_col} columns | "
        f"identity>={args.min_identity:.0%} occupancy>={args.min_occupancy:.0%} "
        f"gaps={'vote' if args.gap_votes else 'ignored'} weighting={args.weighting}\n"
        f"[*] {int(scores.conserved.sum())} conserved columns -> "
        f"{len(runs)} stretches >= {args.min_length} aa -> "
        f"{n_windows} windows of {args.min_length}-{args.max_length} aa",
        file=sys.stderr,
    )

    if args.columns_out:
        with open(args.columns_out, "w") as fh:
            write_columns(fh, scores, ref_pos)

    return 0


if __name__ == "__main__":
    sys.exit(main())
