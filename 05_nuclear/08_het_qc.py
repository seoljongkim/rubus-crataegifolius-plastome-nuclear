#!/usr/bin/env python3
"""08_het_qc.py — Heterozygosity QC of the nuclear SNP set and definition of the strict set (revision R1).

Usage (writes into $RUBUS_BASE/revision_genes/het_qc unless an output dir is given):
  bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\t%MQ0F[\t%GT:%DP:%AD]\n' \
      $RUBUS_BASE/home_wsl_package/nuclear_snp/results/all8.filtered.snps.vcf.gz \
    | python 08_het_qc.py [outdir]


Reads `bcftools query` output on stdin (one line per biallelic SNP):
  CHROM POS QUAL MQ MQ0F  then 8 x  GT:DP:AD
and writes
  ab_hist.tsv        allele-balance histogram of heterozygous calls per sample
  dp_het.tsv         heterozygosity by per-sample depth bin
  mq_het.tsv         heterozygosity by site mapping-quality bin
  strict_gt.tsv.gz   genotype matrix (0/1/2) of sites passing the strict filter
  orig_gt.tsv.gz     genotype matrix of the original polymorphic set (-1 = missing)
  summary.txt        counts
Strict filter: MQ >= 50, per-sample DP 8-40, het allele balance 0.2-0.8, no missing,
polymorphic among the eight cultivars.
"""
import sys, gzip, collections
import os, sys
BASE = os.environ.get("RUBUS_BASE") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

SAMPLES = ["DAHK_2", "DAH_B1", "DB_A1", "DS_A1", "GDCH_1", "GDM_A1", "GDW_5", "HGB_A1"]
N = len(SAMPLES)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "revision_genes", "het_qc")
os.makedirs(OUT, exist_ok=True)
DP_BINS = [(1, 7), (8, 15), (16, 25), (26, 40), (41, 60), (61, 10**9)]
MQ_BINS = [(0, 29), (30, 39), (40, 49), (50, 59), (60, 60)]

def dpbin(d):
    for i, (a, b) in enumerate(DP_BINS):
        if a <= d <= b: return i
    return None
def mqbin(m):
    for i, (a, b) in enumerate(MQ_BINS):
        if a <= m <= b: return i
    return len(MQ_BINS) - 1

ab_hist = [collections.Counter() for _ in range(N)]        # AB bin (0.05) -> count
dp_cnt = [[[0, 0] for _ in DP_BINS] for _ in range(N)]      # [called, het]
mq_cnt = [[[0, 0] for _ in MQ_BINS] for _ in range(N)]
tot = collections.Counter()
strict = gzip.open(f"{OUT}/strict_gt.tsv.gz", "wt")
orig = gzip.open(f"{OUT}/orig_gt.tsv.gz", "wt")

def process(line):
    global tot
    f = line.rstrip("\n").split("\t")
    chrom, pos = f[0], f[1]
    mq = int(float(f[3]))
    gts = []; ok = True; codes = []
    for i in range(N):
        gt, dp, ad = f[5 + i].split(":")
        if gt in ("./.", ".", ".|."):
            codes.append(-1); ok = False; continue
        a, b = gt.replace("|", "/").split("/")
        code = int(a) + int(b)
        codes.append(code)
        dp = int(dp) if dp != "." else 0
        het = code == 1
        bi = dpbin(dp)
        if bi is not None:
            dp_cnt[i][bi][0] += 1; dp_cnt[i][bi][1] += het
        mi = mqbin(mq)
        mq_cnt[i][mi][0] += 1; mq_cnt[i][mi][1] += het
        if het:
            ads = ad.split(",")
            try:
                r, alt = int(ads[0]), int(ads[1])
            except (ValueError, IndexError):
                r, alt = 0, 0
            if r + alt > 0:
                abv = alt / (r + alt)
                ab_hist[i][min(int(abv * 20), 19)] += 1
                if not (0.2 <= abv <= 0.8): ok = False
            else:
                ok = False
        if not (8 <= dp <= 40): ok = False
    tot["sites"] += 1
    called = [c for c in codes if c >= 0]
    poly = len(set(called)) > 1
    if poly:
        tot["orig_poly"] += 1
        orig.write(f"{chrom}\t{pos}\t" + "\t".join(map(str, codes)) + "\n")
    if ok and mq >= 50 and poly:
        tot["strict_poly"] += 1
        strict.write(f"{chrom}\t{pos}\t" + "\t".join(map(str, codes)) + "\n")
    elif ok and mq >= 50:
        tot["strict_mono"] += 1
bad = 0
for line in sys.stdin:
    try:
        process(line)
    except Exception as e:
        bad += 1
        if bad <= 5: sys.stderr.write(f"BAD[{type(e).__name__}: {e}]: {line[:300]}\n")
tot["bad_lines"] = bad
strict.close(); orig.close()

with open(f"{OUT}/ab_hist.tsv", "w") as o:
    o.write("ab_bin_low\t" + "\t".join(SAMPLES) + "\n")
    for k in range(20):
        o.write(f"{k/20:.2f}\t" + "\t".join(str(ab_hist[i][k]) for i in range(N)) + "\n")
with open(f"{OUT}/dp_het.tsv", "w") as o:
    o.write("sample\tdp_bin\tcalled\thet\thet_rate\n")
    for i in range(N):
        for j, (a, b) in enumerate(DP_BINS):
            c, h = dp_cnt[i][j]
            o.write(f"{SAMPLES[i]}\t{a}-{b if b < 10**8 else ''}\t{c}\t{h}\t{h/c if c else float('nan'):.4f}\n")
with open(f"{OUT}/mq_het.tsv", "w") as o:
    o.write("sample\tmq_bin\tcalled\thet\thet_rate\n")
    for i in range(N):
        for j, (a, b) in enumerate(MQ_BINS):
            c, h = mq_cnt[i][j]
            o.write(f"{SAMPLES[i]}\t{a}-{b}\t{c}\t{h}\t{h/c if c else float('nan'):.4f}\n")
with open(f"{OUT}/summary.txt", "w") as o:
    for k, v in tot.items(): o.write(f"{k}\t{v}\n")
print(dict(tot))
