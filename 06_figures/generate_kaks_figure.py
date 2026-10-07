#!/usr/bin/env python3
"""generate_kaks_figure.py — Fig. 3: 플라스톰 코딩서열의 선택압 (리뷰어 R2-4 대응)
(A) 구획별 길이 정규화 치환율  (B) 품종쌍 연결서열 Ka/Ks  (C) 근연종 대비 유전자별 Ka/Ks"""
import csv, json, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

K=f"{BASE}/kaks_output"; OUT=f"{BASE}/figures_revision"
NPG=["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4"]
plt.rcParams.update({"font.family":"serif","font.size":11})
rng=np.random.default_rng(20260819)

# ---- 데이터 ----
part=[]
for l in open(f"{K}/C_coding_noncoding_rates.txt"):
    m=re.match(r"(CDS|tRNA/rRNA|noncoding)\tlength=(\d+)\tvariants=(\d+).*rate_per_kb=([\d.]+)", l)
    if m: part.append((m.group(1), int(m.group(2)), int(m.group(3)), float(m.group(4))))
pairs=[r for r in csv.DictReader(open(f"{K}/A_intraspecific_concat_pairwise.csv")) if r["KaKs"]!="NA"]
# per-pair 95% bootstrap CIs (300 codon-column resamples; values as in Supplementary Table S4)
CI={r["pair"]:(float(r["ci_lo"]),float(r["ci_hi"])) for r in csv.DictReader(open(f"{K}/A_intraspecific_concat_ci.csv"))}
gene={}
for r in csv.DictReader(open(f"{K}/B_interspecific_kaks.csv")):
    if r["KaKs"]!="NA": gene.setdefault(r["species"],[]).append(float(r["KaKs"]))
conc=json.load(open(f"{K}/bootstrap_ci.json"))
SP=["R. takesimensis","R. coreanus","R. chingii"]

fig=plt.figure(figsize=(15.5,5.0))
gs=fig.add_gridspec(1,3,width_ratios=[1.0,1.25,1.15],wspace=0.30,left=.05,right=.985,top=.88,bottom=.20)

# ---- (A) 구획별 치환율 ----
ax=fig.add_subplot(gs[0])
labs=["Protein-coding","tRNA / rRNA","Non-coding"]
rates=[p[3] for p in part]; lens=[p[1] for p in part]; nv=[p[2] for p in part]
b=ax.bar(range(3), rates, 0.6, color=[NPG[3],NPG[2],NPG[0]], edgecolor="white", lw=.8)
for i,(r_,n_,l_) in enumerate(zip(rates,nv,lens)):
    ax.text(i, r_+0.18, f"{r_:.2f}", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text(i, -0.62, f"{n_} variants\n{l_/1000:.1f} kb", ha="center", va="top", fontsize=10, color="#444444")
ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=10.5)
ax.set_ylabel("variants per kb", fontsize=11.5); ax.set_ylim(0, max(rates)*1.20)
ax.annotate("", xy=(0,rates[0]+0.05), xytext=(2,rates[2]+0.05),
            arrowprops=dict(arrowstyle="<->", color="#666666", lw=1.0))
ax.text(1, rates[2]*0.72, "0.23×", ha="center", fontsize=11, color="#444444")
ax.set_title("(a)  Substitution rate by partition", loc="left", fontsize=12.5, fontweight="bold", pad=10)
for s in ("top","right"): ax.spines[s].set_visible(False)

# ---- (B) 품종쌍 Ka/Ks + CI ----
ax=fig.add_subplot(gs[1])
pairs.sort(key=lambda r: float(r["KaKs"]))
y=np.arange(len(pairs)); v=[float(r["KaKs"]) for r in pairs]
lo=[CI.get(r["pair"],(np.nan,np.nan))[0] for r in pairs]; hi=[CI.get(r["pair"],(np.nan,np.nan))[1] for r in pairs]
ax.hlines(y, lo, hi, color="#BBBBBB", lw=1.6, zorder=1)
ax.scatter(v, y, s=26, color=NPG[3], zorder=3)
ax.axvline(1.0, color=NPG[0], ls="--", lw=1.3)
ax.text(0.985, len(pairs)*0.55, "neutrality (Ka/Ks = 1)", color=NPG[0], fontsize=10.5,
        rotation=90, va="center", ha="right")
ax.set_yticks(y); ax.set_yticklabels([r["pair"] for r in pairs], fontsize=10)
ax.set_xlabel("Ka/Ks (concatenated CDS)", fontsize=11.5)
ax.set_xlim(0, 1.12); ax.set_ylim(-1, len(pairs))
ax.set_title("(b)  Among cultivars (22 pairs)", loc="left", fontsize=12.5, fontweight="bold", pad=10)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.grid(axis="x", alpha=.22, ls=":")

# ---- (C) 근연종 대비 유전자별 Ka/Ks ----
ax=fig.add_subplot(gs[2])
for i,sp in enumerate(SP):
    g=np.array(gene[sp]); x=rng.normal(i, 0.075, len(g))
    ax.scatter(x, g, s=17, color=NPG[[1,2,5][i]], alpha=.62, edgecolor="none", zorder=2)
    c=conc[sp]
    ax.errorbar(i, c["point"], yerr=[[c["point"]-c["lo"]],[c["hi"]-c["point"]]], fmt="D",
                ms=7, color="black", ecolor="black", elinewidth=1.5, capsize=5, zorder=4)
    ax.text(i+0.20, c["point"], f"{c['point']:.3f}", fontsize=10.5, va="center")
ax.axhline(1.0, color=NPG[0], ls="--", lw=1.3)
ax.text(2.42, 1.0, "Ka/Ks = 1", color=NPG[0], fontsize=10.5, va="bottom", ha="right")
ax.set_xticks(range(3)); ax.set_xticklabels([f"$\\it{{{s.replace(' ',chr(92)+' ')}}}$\n(n = {len(gene[s])})" for s in SP], fontsize=10)
ax.set_ylabel("Ka/Ks", fontsize=11.5); ax.set_ylim(-0.05, 1.15); ax.set_xlim(-0.5, 2.5)
ax.set_title("(c)  Against congeners, per gene", loc="left", fontsize=12.5, fontweight="bold", pad=10)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", alpha=.22, ls=":")

for ext in ("png","pdf"):
    fig.savefig(f"{OUT}/FigX_selection.{ext}", dpi=300 if ext=="png" else None, bbox_inches="tight")
print(f"저장 완료 · 구획 {len(part)} · 품종쌍 {len(pairs)} · 유전자값 {sum(len(v) for v in gene.values())}")
