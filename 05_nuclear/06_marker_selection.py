#!/usr/bin/env python3
"""
06_marker_selection.py — 품종 식별용 최소 SNP 패널 선정 (RcNu-01 ~ RcNu-05)

선정 기준 (원고 Methods와 동일):
  QUAL >= 200 · 전 품종 DP 8~60 · 결측 없음 · 염색체상 위치
  · ±60 bp 이내에 품종 간 다형 없음 · 121 bp 창 전 구간에서 전 품종 깊이 >= 8x
이후 그리디 집합피복으로 28쌍을 모두 구분하는 최소 조합을 만들고,
가장 가까운 쌍(DAHE-DAHO, DB-GDC)은 이형접합 마커로 보완한다.

입력: cand_raw.tsv  (bcftools query 결과, 아래 CMD 참조)
      all8_poly.bim (plink2 --mac 1 산출물)
      *.md.bam      (samtools depth용)
"""
import collections, itertools, subprocess, json, sys
import numpy as np

CMD = ("bcftools query -i 'QUAL>=100 && MIN(FMT/DP)>=8 && MAX(FMT/DP)<=60 && COUNT(GT=\"mis\")=0' "
       "-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL[\\t%GT]\\n' all8.filtered.snps.vcf.gz > cand_raw.tsv")

LAB   = ["DAHE","DAHO","DB","DS","GDC","GDM","GDW","HGB"]
BAMS  = ["DAHK_2","DAH_B1","DB_A1","DS_A1","GDCH_1","GDM_A1","GDW_5","HGB_A1"]
CODE  = {"0/0":0, "0/1":1, "1/1":2, "1/0":1}
CHROM = {"OZ232399.1","OZ232400.1","OZ232401.1","OZ232402.1","OZ232403.1","OZ232404.1","OZ232405.1"}
BAMDIR, WIN, MINDP = "../bam", 60, 8

def load_poly(bim="all8_poly.bim"):
    d = collections.defaultdict(list)
    for line in open(bim):
        f = line.split(); d[f[0]].append(int(f[3]))
    return {c: np.array(sorted(v)) for c, v in d.items()}

def clear(POLY, c, p, win=WIN):
    a = POLY.get(c)
    if a is None: return False
    i, j = np.searchsorted(a, p-win), np.searchsorted(a, p+win+1)
    return (j - i) <= 1

def cov_ok(c, p, win=WIN, mind=MINDP):
    bams = [f"{BAMDIR}/{s}.md.bam" for s in BAMS]
    out = subprocess.run(["samtools","depth","-a","-Q","20","-r",f"{c}:{p-win}-{p+win}"]+bams,
                         capture_output=True, text=True).stdout.strip().split("\n")
    if len(out) < 2*win+1: return False, 0
    m = 10**9
    for l in out:
        d = [int(x) for x in l.split("\t")[2:]]
        if len(d) < 8: return False, 0
        m = min(m, min(d))
    return m >= mind, m

def read_candidates(POLY, path="cand_raw.tsv", minqual=200, hom_only=True):
    pool = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] not in CHROM or float(f[4]) < minqual: continue
        g = tuple(CODE.get(x, -1) for x in f[5:])
        if -1 in g: continue
        s = set(g)
        if len(s) < 2: continue
        if hom_only and not s <= {0, 2}: continue
        if not clear(POLY, f[0], int(f[1])): continue
        pool.append((f[0], int(f[1]), f[2], f[3], float(f[4]), g))
    return pool

def greedy(cands):
    """28쌍을 최대한 구분하는 최소 조합 (패턴별 후보가 많은 쪽을 우선해 예비 마커 확보)"""
    ALL = list(itertools.combinations(range(8), 2))
    by = collections.defaultdict(list)
    for c in cands: by[c[5]].append(c)
    chosen, cov, avail = [], set(), dict(by)
    while len(cov) < 28 and avail:
        best = max(avail, key=lambda g: (len({p for p in ALL if g[p[0]] != g[p[1]]} - cov), len(avail[g])))
        gain = {p for p in ALL if best[p[0]] != best[p[1]]} - cov
        if not gain: break
        chosen.append(sorted(avail[best], key=lambda y: -y[4])[0]); cov |= gain; del avail[best]
    return chosen, cov

if __name__ == "__main__":
    print("사전 준비:\n  " + CMD + "\n  plink2 --bfile all8 --allow-extra-chr --mac 1 --make-bed --out all8_poly\n")
    POLY = load_poly()
    pool = [c for c in read_candidates(POLY) if cov_ok(c[0], c[1])[0]]
    chosen, cov = greedy(pool)
    print(f"후보 {len(pool)}개 → 최소 패널 {len(chosen)}개로 28쌍 중 {len(cov)}쌍 구분")
    # 남은 쌍은 이형접합 허용 마커로 보완
    ALL = list(itertools.combinations(range(8), 2))
    rest = [p for p in ALL if p not in cov]
    if rest:
        het = read_candidates(POLY, minqual=300, hom_only=False)
        for a, b in rest:
            extra = [c for c in het if c[5][a] != c[5][b] and cov_ok(c[0], c[1])[0]]
            if extra: chosen.append(sorted(extra, key=lambda y: -y[4])[0])
    GT = {0:"AA", 1:"AB", 2:"BB"}
    print("\n=== 최종 패널 ===")
    for i, c in enumerate(chosen, 1):
        print(f"  RcNu-{i:02d} {c[0]}:{c[1]} {c[2]}>{c[3]} QUAL={c[4]:.0f}  " +
              " ".join(f"{l}:{GT[c[5][j]]}" for j, l in enumerate(LAB)))
    codes = ["|".join(GT[c[5][j]] for c in chosen) for j in range(8)]
    print(f"\n바코드 모두 고유: {len(set(codes)) == 8}")
