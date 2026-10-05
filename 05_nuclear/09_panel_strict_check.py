#!/usr/bin/env python3
"""09_panel_strict_check.py — Re-check the marker panel against the strict SNP set and find a replacement for RcNu-04 (revision R1).

Requires strict_gt.tsv.gz from 08_het_qc.py, cand_raw.tsv and all8_poly.bim from 06_marker_selection.py,
and the BAM files; samtools and bcftools must be on PATH. Original motivation: the RcNu-04 of the first
submission (OZ232399.1:125521) rested on a heterozygous call with allele balance 0.12 and failed the strict
filter; OZ232403.1:12029618 was chosen from the ranked list this script prints.

Find a replacement for RcNu-04 that (i) is in the strict SNP set, (ii) has clean allele depths
(het 0.3-0.7, hom with zero minor-allele reads), (iii) separates 'Dongakheuk' from 'Dongakhong',
(iv) satisfies the original window/coverage criteria. Prints ranked candidates."""
import collections, gzip, itertools, subprocess, sys
import os, sys
BASE = os.environ.get("RUBUS_BASE") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np

RES = os.path.join(BASE, "home_wsl_package", "nuclear_snp", "results")
BAMDIR = os.path.join(BASE, "home_wsl_package", "nuclear_snp", "bam")
STRICT = os.path.join(BASE, "revision_genes", "het_qc", "strict_gt.tsv.gz")
LAB = ["DAHE", "DAHO", "DB", "DS", "GDC", "GDM", "GDW", "HGB"]
BAMS = ["DAHK_2", "DAH_B1", "DB_A1", "DS_A1", "GDCH_1", "GDM_A1", "GDW_5", "HGB_A1"]
CODE = {"0/0": 0, "0/1": 1, "1/1": 2, "1/0": 1}
CHROM = {"OZ232399.1", "OZ232400.1", "OZ232401.1", "OZ232402.1", "OZ232403.1", "OZ232404.1", "OZ232405.1"}
WIN, MINDP = 60, 8
PANEL = {("OZ232402.1", 27407814), ("OZ232400.1", 932131), ("OZ232401.1", 5003032), ("OZ232399.1", 1175884)}

poly = collections.defaultdict(list)
for line in open(f"{RES}/all8_poly.bim"):
    f = line.split(); poly[f[0]].append(int(f[3]))
poly = {c: np.array(sorted(v)) for c, v in poly.items()}
def clear(c, p):
    a = poly.get(c)
    if a is None: return False
    i, j = np.searchsorted(a, p - WIN), np.searchsorted(a, p + WIN + 1)
    return (j - i) <= 1

strict = set()
with gzip.open(STRICT, "rt") as f:
    for l in f:
        c, p = l.split("\t", 2)[:2]
        if c in CHROM: strict.add((c, int(p)))
print("strict chromosome sites:", len(strict))

cands = []
for line in open(f"{RES}/cand_raw.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] not in CHROM or float(f[4]) < 300: continue
    key = (f[0], int(f[1]))
    if key not in strict: continue
    g = tuple(CODE.get(x, -1) for x in f[5:])
    if -1 in g or g[0] == g[1]: continue          # must separate DAHE vs DAHO
    if not clear(*key): continue
    cands.append((key, f[2], f[3], float(f[4]), g))
print("candidates separating DAHE/DAHO, strict, clear window:", len(cands))
cands.sort(key=lambda x: -x[3])

# allele depths
regs = os.path.join(os.path.dirname(STRICT), "panel_candidate_regions.txt")
with open(regs, "w") as o:
    for (c, p), *_ in cands[:3000]: o.write(f"{c}\t{p}\n")
out = subprocess.run(["bcftools", "query", "-R", regs, "-f", "%CHROM\t%POS\t%MQ[\t%AD]\n", f"{RES}/all8.filtered.snps.vcf.gz"],
                     capture_output=True, text=True).stdout
ad = {}
for l in out.strip().split("\n"):
    f = l.split("\t"); ad[(f[0], int(f[1]))] = (int(float(f[2])), [tuple(int(x) for x in a.split(",")) for a in f[3:]])
def clean(g, ads):
    for gt, (r, a) in zip(g, ads):
        n = r + a
        if n == 0: return False
        if gt == 1 and not (0.3 <= a / n <= 0.7): return False
        if gt == 0 and a > 0: return False
        if gt == 2 and r > 0: return False
    return True
good = [c for c in cands[:3000] if c[0] in ad and clean(c[4], ad[c[0]][1])]
print("with clean allele depths:", len(good))

def cov_ok(c, p):
    bams = [f"{BAMDIR}/{s}.md.bam" for s in BAMS]
    try:
        o = subprocess.run(["samtools", "depth", "-a", "-Q", "20", "-r", f"{c}:{p-WIN}-{p+WIN}"] + bams,
                           capture_output=True, text=True, timeout=120).stdout.strip().split("\n")
    except subprocess.TimeoutExpired:
        return False, 0
    if len(o) < 2 * WIN + 1: return False, 0
    m = 10 ** 9
    for l in o:
        d = [int(x) for x in l.split("\t")[2:]]
        if len(d) < 8: return False, 0
        m = min(m, min(d))
    return m >= MINDP, m

GT = {0: "AA", 1: "AB", 2: "BB"}
ALL = list(itertools.combinations(range(8), 2))
# how many pairs does the existing 4-marker panel separate with >=2 markers? prefer candidates adding redundancy on DB-GDC / DAHE-DAHO
shown = 0
for key, ref, alt, q, g in good:
    ok, m = cov_ok(*key)
    if not ok: continue
    sep = [p for p in ALL if g[p[0]] != g[p[1]]]
    hets = sum(1 for x in g if x == 1)
    print(f"{key[0]}:{key[1]} {ref}>{alt} QUAL={q:.0f} MQ={ad[key][0]} minDP={m} hets={hets} pairs={len(sep)} "
          + " ".join(f"{l}:{GT[x]}" for l, x in zip(LAB, g)) + "  AD=" + " ".join(f"{r},{a}" for r, a in ad[key][1])
          + ("  [sep DB-GDC]" if g[2] != g[4] else ""))
    shown += 1
    if shown >= 25: break
