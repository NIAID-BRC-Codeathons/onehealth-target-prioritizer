#!/usr/bin/env python3
import importlib.util, io, numpy as np, sys
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

# The worked example from the design review: 6 sequences, 2 stretches at 80%,
# built so that every interesting case appears exactly once.
MOCK = ["MKT-ACDEFGHIKLMNPQ--RSTVWYACDE",   # A  clean
        "MKTLACDEFGHIKLMNPQ--RSTVWYACDE",   # B  clean, different offset
        "M--QACDEFSHIKLMNPQGGRSTVWYACDE",   # C  substitution at col 10
        "MKTWACDEFGHIKLMNPQ--RSTVWYASDE",   # D  substitution at col 28
        "MKT-ACDEFGHIKLMNPQ--RST---ACDE",   # E  internal deletion
        "MKTRACDEFGHIKLMNPQ------------"]   # F  truncated
mc = codes_from(MOCK)
ms = score(MOCK, min_identity=0.80)
mruns = cr.maximal_runs(ms.conserved, 10)
mrp = cr.reference_map(mc)

print("\n[9] per-sequence identity within a stretch")
check("two stretches at 80%", mruns, [(4, 17), (20, 29)])

def ident(i, s, e, gap_votes=False):
    v = cr.subsequence_identity(mc[i, s:e+1], ms.top_code[s:e+1],
                                ms.has_aa[s:e+1], gap_votes)
    return None if v is None else round(v, 4)

check("exact match", ident(0, 4, 17), 1.0)
check("one substitution -> 13/14", ident(2, 4, 17), 0.9286)
check("deletion, gaps ignored -> 7/7", ident(4, 20, 29), 1.0)
check("deletion, gaps vote    -> 7/10", ident(4, 20, 29, True), 0.7)
check("all gaps -> None", ident(5, 20, 29), None)

ac = codes_from(["AAAA", "AAAA", "AAXA"])
asc = score(["AAAA", "AAAA", "AAXA"], min_identity=0.6)
check("ambiguity code scores as a mismatch, never as the consensus",
      round(cr.subsequence_identity(ac[2], asc.top_code, asc.has_aa, False), 4), 0.75)

check("span_label as one cell", cr.span_label(mrp[0], 20, 29), "18-27")
check("span_label when absent", cr.span_label(mrp[5], 20, 29), "NA")
check("reference_map vectorises over rows",
      list(map(int, mrp[2])), list(map(int, cr.reference_map(mc[2]))))

print("\n[10] matrix and fasta writers")
maln = cr.Alignment(
    ids=[f"sp_{c}" for c in "ABCDEF"], codes=mc,
    chars=np.array([bytearray(s, "ascii") for s in MOCK], dtype=np.uint8))

buf = io.StringIO()
cr.write_matrix(buf, maln, mruns, ms, mrp, False)
rows = [r.split("\t") for r in buf.getvalue().rstrip("\n").split("\n")]
check("header 1 repeats each span 3x",
      rows[0], ["region", "5-18", "5-18", "5-18", "21-30", "21-30", "21-30"])
check("header 2 names the fields", rows[1][1:4], ["aligned", "residues", "identity"])
check("consensus row is present", rows[2][:2], ["consensus", "ACDEFGHIKLMNPQ"])
check("2 headers + consensus + 6 sequences", len(rows), 9)
check("aligned cell keeps alignment gaps", rows[7][4], "RST---ACDE")
check("residue cell uses the sequence's own numbering", rows[7][5], "18-24")
check("all-gap stretch -> NA, not a fake span", rows[8][5:7], ["NA", "NA"])

buf = io.StringIO()
n_rec, n_skip = cr.write_fasta(buf, maln, mruns, ms, mrp, False)
fa = buf.getvalue().rstrip("\n").split("\n")
check("records written", n_rec, 11)
check("all-gap sequence skipped", n_skip, 1)
check("2 lines per record", len(fa), 22)
check("grouped by stretch", fa[0].split()[0], ">sp_A_S1")
check("header carries region, residues and identity",
      fa[4], ">sp_C_S1 region=5-18 residues=3-16 identity=0.9286")
check("gaps stripped, so 7 aa not 10", fa[-1], "RSTACDE")

print("\n" + ("ALL PASS" if ok else "FAILURES PRESENT"))
sys.exit(0 if ok else 1)
