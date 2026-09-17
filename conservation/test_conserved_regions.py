#!/usr/bin/env python3
import importlib.util, numpy as np, sys
spec = importlib.util.spec_from_file_location("cr", "conserved_regions.py")
cr = importlib.util.module_from_spec(spec); sys.modules["cr"] = cr; spec.loader.exec_module(cr)

ok = True
def check(label, got, want):
    global ok
    good = got == want
    ok &= good
    print(f"  {'PASS' if good else 'FAIL'}  {label}: got {got!r}, want {want!r}")

def codes_from(seqs):
    raw = np.array([bytearray(s, "ascii") for s in seqs], dtype=np.uint8)
    return cr._CODE_OF[raw]

def score(seqs, **kw):
    kw.setdefault("min_identity", 0.95); kw.setdefault("min_occupancy", 0.5)
    kw.setdefault("gap_votes", False)
    c = codes_from(seqs)
    return cr.score_columns(c, np.ones(len(seqs)), **kw)

print("\n[1] tiling math: a run of exactly 40 columns, windows 10-15")
check("window count", len(list(cr.tile(0, 39, 10, 15))), 31+30+29+28+27+26)
check("run of 10 -> single 10mer", list(cr.tile(0, 9, 10, 15)), [(0, 9)])
check("run of 9  -> nothing", list(cr.tile(0, 8, 10, 15)), [])
check("first window", list(cr.tile(0, 39, 10, 15))[0], (0, 9))
check("last window", list(cr.tile(0, 39, 10, 15))[-1], (30, 39))
check("all windows within bounds",
      all(0 <= s and e <= 39 and 10 <= e-s+1 <= 15 for s, e in cr.tile(0,39,10,15)), True)

print("\n[2] gap voting: column of 9 A + 1 gap (occupancy 0.9)")
col = ["A"]*9 + ["-"]
off = score(col, gap_votes=False); on = score(col, gap_votes=True)
check("gaps ignored -> identity", round(float(off.identity[0]), 4), 1.0)
check("gaps ignored -> conserved", bool(off.conserved[0]), True)
check("gaps vote    -> identity", round(float(on.identity[0]), 4), 0.9)
check("gaps vote    -> conserved at 95%", bool(on.conserved[0]), False)

print("\n[3] occupancy guard: 3 of 20 sequences have A, rest gaps")
col = ["A"]*3 + ["-"]*17
s = score(col, gap_votes=False)
check("identity is 100% among residues", round(float(s.identity[0]), 4), 1.0)
check("occupancy", round(float(s.occupancy[0]), 4), 0.15)
check("guard blocks it", bool(s.conserved[0]), False)
check("no guard -> would pass", bool(score(col, min_occupancy=0.0).conserved[0]), True)

print("\n[4] ambiguity codes occupy but never win: 6 X + 4 A")
s = score(["X"]*6 + ["A"]*4)
check("consensus is A not X", s.consensus[0], "A")
check("X counts in denominator -> 4/10", round(float(s.identity[0]), 4), 0.4)
check("occupancy counts X as occupied", round(float(s.occupancy[0]), 4), 1.0)

print("\n[5] degenerate columns")
check("all gaps -> consensus '-'", score(["-"]*10).consensus[0], "-")
check("all gaps -> identity 0", float(score(["-"]*10).identity[0]), 0.0)
check("all gaps -> not conserved", bool(score(["-"]*10).conserved[0]), False)
check("all X    -> consensus '-'", score(["X"]*10).consensus[0], "-")
check("lowercase treated as residue", score(["a"]*10).consensus[0], "A")

print("\n[6] threshold is inclusive, and 95 == 0.95")
check("exactly 95% passes", bool(score(["A"]*19 + ["C"], gap_votes=True).conserved[0]), True)
check("94.7% fails", bool(score(["A"]*18 + ["C"], gap_votes=True).conserved[0]), False)
check("fraction('95')", cr.fraction("95"), 0.95)
check("fraction('0.95')", cr.fraction("0.95"), 0.95)

print("\n[7] maximal_runs is strictly contiguous (no bridging)")
m = np.array([1,1,1,0,1,1,1], dtype=bool)
check("split at the gap", cr.maximal_runs(m, 1), [(0,2),(4,6)])
check("min_length filter", cr.maximal_runs(m, 4), [])

print("\n[8] reference coordinate mapping")
rp = cr.reference_map(codes_from(["--ACD-EF"])[0])
check("positions", list(map(int, rp)), [0,0,1,2,3,0,4,5])
check("span skips leading gaps", cr.reference_span(rp, 0, 3), ("1","2"))
check("span cols 0-4 = residues 1-3", cr.reference_span(rp, 0, 4), ("1","3"))
check("span across internal gap", cr.reference_span(rp, 4, 6), ("3","4"))
check("all-gap span", cr.reference_span(rp, 0, 1), ("NA","NA"))

print("\n" + ("ALL PASS" if ok else "FAILURES PRESENT"))
sys.exit(0 if ok else 1)
