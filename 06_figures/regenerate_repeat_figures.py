#!/usr/bin/env python3
"""regenerate_repeat_figures.py — Fig 2C/2D/2E 재생성.
원본 패널 스타일(NPG 팔레트, serif, 범례 하단) 유지.
2026-09-14 수정: (1) 품종명 'Dongakhong' 오타 교정,
(2) IR 유래 회문 반복(IRa/IRb 쌍, 25,781~25,786 bp)을 분산반복 집계에서 제외 —
Fig 2E의 길이 범위 밖이어서 Table 2 합계와 어긋나던 문제를 해소한다."""
import os, csv, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = os.environ.get("RUBUS_BASE", "/home/spenkim/rubus_cp")
TRF  = f"{BASE}/trf_output"
REP  = f"{BASE}/reputer_output"
OUT  = f"{BASE}/figures_revision"
os.makedirs(OUT, exist_ok=True)

NPG   = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]
NAME  = {"DB":"Dabok","DS":"Daesung","DAHE":"Dongakheuk","DAHO":"Dongakhong",
         "GDC":"Geumdongcheon","GDW":"Geumdongwang","GDM":"Goldmoon","HGB":"Hwangeumboll"}
ORDER = ["DB","DS","DAHE","DAHO","GDC","GDW","GDM","HGB"]
plt.rcParams.update({"font.family":"serif","font.size":13,
                     "axes.spines.top":False,"axes.spines.right":False})

# ---------------- 데이터 ----------------
def tandem(v):
    f = glob.glob(f"{TRF}/{v}_std.fasta.*.dat")[0]
    L = []
    for line in open(f):
        p = line.split()
        if len(p) >= 15 and p[0].isdigit() and p[1].isdigit():
            L.append(int(p[1]) - int(p[0]) + 1)
    return L

def dispersed(v):
    """IR 유래 반복은 제외한다."""
    rows = list(csv.DictReader(open(f"{REP}/{v}_dispersed_repeats.tsv"), delimiter="\t"))
    return [r for r in rows if r["IR_derived"].strip().lower() not in ("yes","y","true")]

TAND = {v: tandem(v) for v in ORDER}
DISP = {v: dispersed(v) for v in ORDER}

# ---------------- Fig 2C: 탠덤반복 길이 구간 누적 ----------------
TCLS = [("20–39",20,39), ("40–59",40,59), ("60–89",60,89), ("≥90",90,10**9)]
fig, ax = plt.subplots(figsize=(8.8, 6.3))
x = np.arange(len(ORDER)); bottom = np.zeros(len(ORDER))
for i,(lab,lo,hi) in enumerate(TCLS):
    vals = np.array([sum(1 for L in TAND[v] if lo <= L <= hi) for v in ORDER])
    ax.bar(x, vals, 0.62, bottom=bottom, color=NPG[i], edgecolor="white", linewidth=0.6, label=lab)
    bottom += vals
for i, tot in enumerate(bottom):
    ax.text(i, tot + 0.7, str(int(tot)), ha="center", va="bottom", fontsize=12, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([NAME[v] for v in ORDER], rotation=45, ha="right")
ax.set_ylabel("Number of tandem repeats", fontsize=15)
ax.set_ylim(0, max(bottom)*1.18)
ax.legend(handles=[Patch(facecolor=NPG[i], label=l) for i,(l,_,_) in enumerate(TCLS)],
          loc="upper center", bbox_to_anchor=(0.5,-0.28), ncol=4, fontsize=12, frameon=False)
plt.tight_layout()
fig.savefig(f"{OUT}/Fig2C_tandem_length.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2C_tandem_length.pdf", bbox_inches="tight"); plt.close(fig)

# ---------------- Fig 2D: 분산반복 유형 (IR 제외) ----------------
TYPES = [("Forward","F"), ("Palindromic","P"), ("Reverse","R")]
fig, ax = plt.subplots(figsize=(8.8, 6.3))
w = 0.26
for i,(lab,code) in enumerate(TYPES):
    vals = [sum(1 for r in DISP[v] if r["Type"] == code) for v in ORDER]
    ax.bar(x + (i-1)*w, vals, w, color=NPG[i], edgecolor="white", linewidth=0.6, label=lab)
for i, v in enumerate(ORDER):
    ax.text(i, len(DISP[v])*0 + max(sum(1 for r in DISP[v] if r["Type"]==c) for _,c in TYPES) + 0.35,
            str(len(DISP[v])), ha="center", va="bottom", fontsize=12, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([NAME[v] for v in ORDER], rotation=45, ha="right")
ax.set_ylabel("Number of repeats", fontsize=15)
ax.set_ylim(0, 14)
ax.legend(handles=[Patch(facecolor=NPG[i], label=l) for i,(l,_) in enumerate(TYPES)],
          loc="upper center", bbox_to_anchor=(0.5,-0.28), ncol=3, fontsize=12, frameon=False)
plt.tight_layout()
fig.savefig(f"{OUT}/Fig2D_dispersed_type.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2D_dispersed_type.pdf", bbox_inches="tight"); plt.close(fig)

# ---------------- Fig 2E: 분산반복 길이 분포 (IR 제외) ----------------
DCLS = [("30–34",30,34), ("35–39",35,39), ("40–49",40,49), ("50–62",50,62)]
alllen = [int(r["Length"]) for v in ORDER for r in DISP[v]]
fig, ax = plt.subplots(figsize=(7.0, 6.3))
vals = [sum(1 for L in alllen if lo <= L <= hi) for _,lo,hi in DCLS]
ax.bar(range(len(DCLS)), vals, 0.55, color=NPG[1], edgecolor="white", linewidth=0.6)
for i,c in enumerate(vals):
    ax.text(i, c + max(vals)*0.015, str(c), ha="center", va="bottom", fontsize=12, fontweight="bold")
ax.set_xticks(range(len(DCLS))); ax.set_xticklabels([f"{l} bp" for l,_,_ in DCLS])
ax.set_xlabel("Repeat length (bp)", fontsize=15)
ax.set_ylabel("Number of repeats", fontsize=15)
ax.set_ylim(0, max(vals)*1.14)
plt.tight_layout()
fig.savefig(f"{OUT}/Fig2E_dispersed_length.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2E_dispersed_length.pdf", bbox_inches="tight"); plt.close(fig)

# ---------------- 요약 ----------------
print(f"{'cultivar':16s} {'tandem':>7s} {'dispersed':>10s}  (F/P/R)")
for v in ORDER:
    c = {t: sum(1 for r in DISP[v] if r["Type"]==t) for t in "FPR"}
    print(f"{NAME[v]:16s} {len(TAND[v]):7d} {len(DISP[v]):10d}  ({c['F']}/{c['P']}/{c['R']})")
print(f"\n탠덤 합계 {sum(len(TAND[v]) for v in ORDER)} | 분산 합계(IR 제외) {sum(len(DISP[v]) for v in ORDER)}")
print(f"Fig 2E 구간 합계 {sum(vals)}  ← 분산 합계와 일치해야 한다")
print(f"출력: {OUT}/Fig2C_tandem_length, Fig2D_dispersed_type, Fig2E_dispersed_length (.png/.pdf)")
