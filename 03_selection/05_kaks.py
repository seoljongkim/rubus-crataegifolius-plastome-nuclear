#!/usr/bin/env python3
"""
05_kaks.py — 선택압 분석 (HEB 리뷰어 R2-4 대응)  [v2]

리뷰어 지적: "intergenic 변이 83%를 'strong purifying selection'의 근거로 해석했으나
              정규화된 coding/noncoding 치환율이나 syn/nonsyn 분석을 수행하지 않았다"

리뷰어가 명시한 두 방법을 모두 수행:
  Part A. 종내(8품종) — 유전자별 동의/비동의 치환 부위 수 + 농축 CDS 연결서열 Ka/Ks
  Part B. 종간 — R. takesimensis / R. coreanus / R. chingii 대비 유전자별·전체 Ka/Ks
          (종내 분기가 극히 얕아 Ka/Ks 추정이 불안정하므로, 안정적인 선택압 값은 종간에서 얻음)
  Part C. coding vs noncoding 치환율 정규화 (길이 보정)

v2 수정: (1) cal_dn_ds가 종결코돈에서 KeyError → 말단 종결코돈 제거 후 계산
         (2) InDel의 Details 형식(gap_in=)이 SNP와 달라 Part C에서 누락되던 것 처리
         (3) Part A를 쌍 누적 계수 대신 '변이 부위 수' 기준으로 변경 (그룹 크기 편향 제거)
"""
import csv, itertools, warnings, os
from collections import defaultdict, Counter
warnings.filterwarnings("ignore")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

OUT  = f"{BASE}/kaks_output"
VARS = ["DAHE","DAHO","DB","DS","GDC","GDM","GDW","HGB"]      # GB048/GB066 제외
EXCLUDE = {"GB048","GB066"}
STOPS = {"TAA","TAG","TGA"}
OUTGROUPS = [("R. takesimensis","NC_037991.1"), ("R. coreanus","MT478114"), ("R. chingii","NC_050296.1")]

# ---------------------------------------------------------------- helpers
def load_cds(gbfile):
    rec = SeqIO.read(gbfile, "genbank")
    out = {}
    for f in rec.features:
        if f.type != "CDS": continue
        g = f.qualifiers.get("gene", f.qualifiers.get("locus_tag", ["?"]))[0]
        if g in out: continue                              # IR 중복 사본 스킵
        s = str(f.extract(rec.seq)).upper()
        if len(s) % 3 or len(s) < 60: continue
        out[g] = s
    return out, len(rec.seq), rec

def trim_stop(s):
    """말단 종결코돈 제거 — cal_dn_ds가 종결코돈에서 KeyError를 냄."""
    return s[:-3] if len(s) >= 6 and s[-3:] in STOPS else s

def clean_pair(s1, s2):
    """길이 맞추고, 갭/N/종결코돈 포함 코돈 컬럼 제거 (pairwise deletion)."""
    s1, s2 = trim_stop(s1), trim_stop(s2)
    n = min(len(s1), len(s2)) // 3 * 3
    o1 = o2 = ""
    for i in range(0, n, 3):
        c1, c2 = s1[i:i+3], s2[i:i+3]
        if any(ch not in "ACGT" for ch in c1 + c2): continue
        if c1 in STOPS or c2 in STOPS: continue
        o1 += c1; o2 += c2
    return o1, o2

def count_syn_nonsyn(s1, s2):
    syn = nonsyn = 0
    for i in range(0, len(s1), 3):
        c1, c2 = s1[i:i+3], s2[i:i+3]
        if c1 == c2: continue
        d = sum(1 for x, y in zip(c1, c2) if x != y)
        if str(Seq(c1).translate(table=11)) == str(Seq(c2).translate(table=11)): syn += d
        else: nonsyn += d
    return syn, nonsyn

def kaks(s1, s2, method="NG86"):
    try:
        dn, ds = cal_dn_ds(CodonSeq(s1), CodonSeq(s2), method=method)
        return dn, ds
    except Exception:
        return None, None

def fmt(x): return f"{x:.6f}" if isinstance(x, float) else "NA"
def ratio(dn, ds): return f"{dn/ds:.4f}" if (dn is not None and ds not in (None, 0) and ds > 0) else "NA"

# ---------------------------------------------------------------- Part A
print("=" * 80)
print("Part A. 종내 선택압 — 8품종 CDS의 동의/비동의 변이 부위")
print("=" * 80)

cds = {v: load_cds(f"{BASE}/bankit_gb/{v}.gb")[0] for v in VARS}
genes = sorted(set.intersection(*(set(cds[v]) for v in VARS)))
eq_len = [g for g in genes if len({len(cds[v][g]) for v in VARS}) == 1]
print(f"8품종 공통 CDS {len(genes)}개 (길이 동일 {len(eq_len)}개 — 인델 보유 {len(genes)-len(eq_len)}개는 부위 계수에서 제외)\n")

gene_rows = []
for g in eq_len:
    L = len(cds[VARS[0]][g]); syn_sites = nonsyn_sites = 0
    for i in range(0, L - 2, 3):
        codons = {cds[v][g][i:i+3] for v in VARS}
        if len(codons) == 1: continue
        if any(c in STOPS or any(ch not in "ACGT" for ch in c) for c in codons): continue
        aas = {str(Seq(c).translate(table=11)) for c in codons}
        var_pos = sum(1 for p in range(3) if len({c[p] for c in codons}) > 1)
        if len(aas) == 1: syn_sites += var_pos
        else: nonsyn_sites += var_pos
    if syn_sites or nonsyn_sites:
        gene_rows.append(dict(gene=g, len_bp=L, syn_sites=syn_sites, nonsyn_sites=nonsyn_sites,
                              N_over_S=(f"{nonsyn_sites/syn_sites:.2f}" if syn_sites else "inf(S=0)")))

with open(f"{OUT}/A_intraspecific_gene_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["gene","len_bp","syn_sites","nonsyn_sites","N_over_S"]); w.writeheader(); w.writerows(gene_rows)

S = sum(r["syn_sites"] for r in gene_rows); N = sum(r["nonsyn_sites"] for r in gene_rows)
print(f"변이 보유 CDS: {len(gene_rows)}개 / 동일 CDS: {len(eq_len)-len(gene_rows)}개")
print(f"{'gene':10s} {'bp':>6s} {'syn':>5s} {'nonsyn':>7s} {'N/S':>10s}")
for r in sorted(gene_rows, key=lambda x: -(x["syn_sites"] + x["nonsyn_sites"]))[:12]:
    print(f"{r['gene']:10s} {r['len_bp']:6d} {r['syn_sites']:5d} {r['nonsyn_sites']:7d} {r['N_over_S']:>10s}")
print(f"\n합계: 동의 부위 {S}, 비동의 부위 {N}, N/S = {N/S:.2f}" if S else f"\n합계: 동의 {S}, 비동의 {N}")

# 농축 CDS 연결서열 pairwise Ka/Ks (전체 28쌍)
print("\n--- 전체 CDS 연결서열 pairwise Ka/Ks (28쌍) ---")
pair_rows = []
for v1, v2 in itertools.combinations(VARS, 2):
    c1 = c2 = ""
    for g in eq_len:
        a, b = clean_pair(cds[v1][g], cds[v2][g]); c1 += a; c2 += b
    s, n = count_syn_nonsyn(c1, c2)
    dn, ds = kaks(c1, c2)
    pair_rows.append(dict(pair=f"{v1}-{v2}", aligned_bp=len(c1), syn=s, nonsyn=n,
                          Ka=fmt(dn), Ks=fmt(ds), KaKs=ratio(dn, ds)))
with open(f"{OUT}/A_intraspecific_concat_pairwise.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(pair_rows[0].keys())); w.writeheader(); w.writerows(pair_rows)
vals = [float(r["KaKs"]) for r in pair_rows if r["KaKs"] != "NA"]
print(f"Ka/Ks 산출 {len(vals)}/28쌍" + (f" | 범위 {min(vals):.4f}~{max(vals):.4f} | 중앙값 {sorted(vals)[len(vals)//2]:.4f}" if vals else ""))
for r in sorted(pair_rows, key=lambda x: (x["KaKs"] == "NA", x["pair"]))[:6]:
    print(f"  {r['pair']:12s} syn={r['syn']:3d} nonsyn={r['nonsyn']:3d}  Ka={r['Ka']:>8s} Ks={r['Ks']:>8s} Ka/Ks={r['KaKs']:>7s}")

# ---------------------------------------------------------------- Part B
print("\n" + "=" * 80)
print("Part B. 종간 선택압 — R. crataegifolius(DAHO) vs 근연종")
print("=" * 80)

all_b = []
for sp, acc in OUTGROUPS:
    other = load_cds(f"{OUT}/{acc}.gb")[0]
    shared = sorted(set(cds["DAHO"]) & set(other))
    rows, cc1, cc2 = [], "", ""
    for g in shared:
        if len(cds["DAHO"][g]) != len(other[g]): continue          # 인델 보유 유전자는 제외
        a, b = clean_pair(cds["DAHO"][g], other[g])
        if len(a) < 90: continue
        cc1 += a; cc2 += b
        s, n = count_syn_nonsyn(a, b)
        dn, ds = kaks(a, b)
        rows.append(dict(species=sp, gene=g, aligned_bp=len(a), syn=s, nonsyn=n,
                         Ka=fmt(dn), Ks=fmt(ds), KaKs=ratio(dn, ds)))
    dn, ds = kaks(cc1, cc2)
    s, n = count_syn_nonsyn(cc1, cc2)
    v = [float(r["KaKs"]) for r in rows if r["KaKs"] != "NA"]
    print(f"\n[{sp} / {acc}] 비교 유전자 {len(rows)}개, 연결서열 {len(cc1):,} bp")
    print(f"  전체 연결서열: 동의 {s}, 비동의 {n} | Ka={fmt(dn)} Ks={fmt(ds)} **Ka/Ks={ratio(dn,ds)}**")
    if v:
        print(f"  유전자별 Ka/Ks 산출 {len(v)}개 | 중앙값 {sorted(v)[len(v)//2]:.4f} | "
              f"<1 정화선택 {sum(1 for x in v if x < 1)}개 | >1 {sum(1 for x in v if x > 1)}개 "
              f"{[r['gene'] for r in rows if r['KaKs'] != 'NA' and float(r['KaKs']) > 1]}")
    all_b += rows

with open(f"{OUT}/B_interspecific_kaks.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(all_b[0].keys())); w.writeheader(); w.writerows(all_b)

# ---------------------------------------------------------------- Part C
print("\n" + "=" * 80)
print("Part C. coding vs noncoding 치환율 정규화 (길이 보정)")
print("=" * 80)

_, glen, rec = load_cds(f"{BASE}/bankit_gb/DAHO.gb")
cds_pos, rna_pos = set(), set()
for f in rec.features:
    if f.type == "CDS":
        for p in f.location.parts: cds_pos.update(range(int(p.start) + 1, int(p.end) + 1))
for f in rec.features:
    if f.type in ("tRNA", "rRNA"):
        for p in f.location.parts: rna_pos.update(range(int(p.start) + 1, int(p.end) + 1))
rna_pos -= cds_pos
non_len = glen - len(cds_pos) - len(rna_pos)

def polymorphic_in8(details, vtype):
    """SNP: 'T:DAHE,...; A:DB,...'  /  InDel: 'gap_in=DB,DS,...; seq=...'"""
    if vtype == "InDel" or details.startswith("gap_in="):
        blk = details.split(";")[0]
        gap = {v.strip() for v in blk.split("=", 1)[1].split(",") if v.strip()} - EXCLUDE
        return 0 < len(gap) < len(VARS)                    # 일부만 갭이면 8품종 내 다형
    alleles = []
    for b in details.split(";"):
        if ":" not in b: continue
        al, vs = b.split(":", 1)
        if [v.strip() for v in vs.split(",") if v.strip() and v.strip() not in EXCLUDE]:
            alleles.append(al.strip())
    return len(alleles) >= 2

cnt, tot = Counter(), Counter()
for r in csv.DictReader(open(f"{BASE}/variant_output/all_variants.csv")):
    if not polymorphic_in8(r["Details"], r["Type"]): continue
    tot[r["Type"]] += 1
    p = int(r["Ref_pos"])
    reg = "CDS" if p in cds_pos else ("tRNA/rRNA" if p in rna_pos else "noncoding")
    cnt[reg] += 1; cnt[(reg, r["Type"])] += 1

T = sum(tot.values())
print(f"8품종 기준 변이 {T}개 (SNP {tot['SNP']}, InDel {tot['InDel']}, Complex {tot['Complex']})\n")
base = cnt["noncoding"] / (non_len / 1000)
print(f"{'구획':12s} {'길이(bp)':>10s} {'게놈비':>7s} {'변이':>6s} {'SNP':>5s} {'InDel':>6s} {'변이비':>7s} {'변이/kb':>9s} {'상대율':>8s}")
for lab, ln in (("CDS", len(cds_pos)), ("tRNA/rRNA", len(rna_pos)), ("noncoding", non_len)):
    rate = cnt[lab] / (ln / 1000)
    print(f"{lab:12s} {ln:10,d} {100*ln/glen:6.1f}% {cnt[lab]:6d} {cnt[(lab,'SNP')]:5d} {cnt[(lab,'InDel')]:6d} "
          f"{100*cnt[lab]/T:6.1f}% {rate:9.3f} {rate/base:7.2f}x")
print(f"\n→ CDS 치환율은 noncoding의 {(cnt['CDS']/(len(cds_pos)/1000))/base:.2f}배 "
      f"(길이 보정 후). tRNA/rRNA는 {(cnt['tRNA/rRNA']/(len(rna_pos)/1000))/base:.2f}배.")

with open(f"{OUT}/C_coding_noncoding_rates.txt", "w") as fh:
    fh.write(f"Genome (DAHO) length\t{glen}\nVariants among 8 cultivars\t{T} (SNP {tot['SNP']}, InDel {tot['InDel']}, Complex {tot['Complex']})\n")
    for lab, ln in (("CDS", len(cds_pos)), ("tRNA/rRNA", len(rna_pos)), ("noncoding", non_len)):
        fh.write(f"{lab}\tlength={ln}\tvariants={cnt[lab]}\tSNP={cnt[(lab,'SNP')]}\tInDel={cnt[(lab,'InDel')]}\t"
                 f"rate_per_kb={cnt[lab]/(ln/1000):.4f}\trelative={((cnt[lab]/(ln/1000))/base):.4f}\n")

print(f"\n결과 파일: {OUT}/A_*.csv, B_interspecific_kaks.csv, C_coding_noncoding_rates.txt")
