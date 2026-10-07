#!/usr/bin/env python3
"""regenerate_ssr_figures.py — 보수적 MISA threshold(5-3, 6-3) 결과로 Fig 2A/2B 재생성.
원본 패널 스타일(NPG 팔레트, serif, 범례 별도) 유지. 원본은 보존하고 figures_revision/ 에 출력."""
import os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

SRC = f"{BASE}/home_wsl_package/misa_new"
OUT = f"{BASE}/figures_revision"
os.makedirs(OUT, exist_ok=True)

NPG   = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]
LABEL = ["Mono", "Di", "Tri", "Tetra", "Penta", "Hexa"]
NAME  = {"DB":"Dabok","DS":"Daesung","DAHE":"Dongakheuk","DAHO":"Dongakhong",
         "GDC":"Geumdongcheon","GDW":"Geumdongwang","GDM":"Goldmoon","HGB":"Hwangeumboll"}
ORDER = ["DB","DS","DAHE","DAHO","GDC","GDW","GDM","HGB"]     # 원본 그림과 동일한 알파벳 순

plt.rcParams.update({"font.family":"serif", "font.size":22,
                     "axes.spines.top":False, "axes.spines.right":False})

def parse(v):
    s = open(f"{SRC}/{v}_std.fasta.statistics").read()
    blk = re.search(r"Unit size\tNumber of SSRs\n(.*?)\n\n", s, re.S).group(1)
    d = {int(a): int(b) for a, b in (l.split("\t") for l in blk.strip().split("\n"))}
    return [d.get(u, 0) for u in range(1, 7)]

counts = {v: parse(v) for v in ORDER}

# ---------------- Fig 2A: 누적 막대 ----------------
fig, ax = plt.subplots(figsize=(8.8, 6.3))
x = np.arange(len(ORDER)); bottom = np.zeros(len(ORDER))
for i in range(6):
    vals = np.array([counts[v][i] for v in ORDER])
    ax.bar(x, vals, 0.62, bottom=bottom, color=NPG[i], edgecolor="white", linewidth=0.6)
    bottom += vals
for i, tot in enumerate(bottom):
    ax.text(i, tot + 1.2, str(int(tot)), ha="center", va="bottom", fontsize=20, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([NAME[v] for v in ORDER], rotation=45, ha="right")
ax.set_ylabel("Number of SSRs", fontsize=24)
ax.set_ylim(0, max(bottom) * 1.16)
plt.tight_layout()
fig.savefig(f"{OUT}/Fig2A_SSR_type.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2A_SSR_type.pdf", bbox_inches="tight"); plt.close(fig)

# ---------------- Fig 2B: 평균 조성 파이 ----------------
avg = np.array([np.mean([counts[v][i] for v in ORDER]) for i in range(6)])
pct = 100 * avg / avg.sum()
keep = [i for i in range(6) if avg[i] > 0]          # hexa는 0건이라 파이에서 제외
fig, ax = plt.subplots(figsize=(7.4, 6.3))
w, _ = ax.pie(avg[keep], colors=[NPG[i] for i in keep], startangle=90, counterclock=False,
              wedgeprops={"edgecolor":"white","linewidth":1.5}, radius=1.0)
ax.axis("equal")
# labels outside the pie with leader lines (no separate legend needed)
import math
for wd, i in zip(w, keep):
    ang = math.radians((wd.theta1 + wd.theta2) / 2)
    if LABEL[i] == "Penta": ang = math.radians(97)      # nudge the 0.7 % wedge label clear of 'Tetra'
    x, y = math.cos(ang), math.sin(ang)
    ha = "left" if x >= 0 else "right"
    ax.annotate(f"{LABEL[i]} {pct[i]:.1f}%", xy=(0.96*x, 0.96*y), xytext=(1.32*x, 1.28*y),
                ha=ha, va="center", fontsize=20,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=1.0,
                                connectionstyle=f"angle,angleA=0,angleB={math.degrees(ang):.0f}"))
ax.text(0, -1.42, f"Hexa 0.0% (none detected)", ha="center", va="top", fontsize=17, color="#555555")
ax.set_xlim(-1.9, 1.9); ax.set_ylim(-1.55, 1.45)
plt.tight_layout()
fig.savefig(f"{OUT}/Fig2B_SSR_proportion.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2B_SSR_proportion.pdf", bbox_inches="tight"); plt.close(fig)

print(f"{'cultivar':16s} " + "".join(f"{l:>7s}" for l in LABEL) + f"{'total':>8s}")
for v in ORDER:
    print(f"{NAME[v]:16s} " + "".join(f"{c:7d}" for c in counts[v]) + f"{sum(counts[v]):8d}")
print("\n평균 조성(%): " + ", ".join(f"{LABEL[i]} {pct[i]:.1f}" for i in range(6)))
print(f"출력: {OUT}/Fig2A_SSR_type.png(.pdf), Fig2B_SSR_proportion.png(.pdf)")

# ---------------- pptx 조립용 분리 이미지 ----------------
# Figure2.pptx는 Fig2B를 '파이(정사각 crop)'와 '범례 띠(하단 crop)' 두 도형으로 배치한다.
# 새 이미지는 내부 비율이 달라 기존 crop이 맞지 않으므로, 두 요소를 각각 별도 이미지로 출력하고
# pptx 쪽에서는 crop을 제거한다. (그림 12 = 1.73x1.73in 정사각, 그림 31 = 1.83x0.30in ≒ 6.02:1)

# 파이 단독 (정사각)
fig, ax = plt.subplots(figsize=(6.0, 6.0))
ax.pie(avg[keep], colors=[NPG[i] for i in keep], startangle=90,
       autopct=lambda p: f"{p:.1f}%" if p >= 3 else "", pctdistance=0.72,
       wedgeprops={"edgecolor":"white","linewidth":1.2}, textprops={"fontsize":15})
ax.axis("equal")
plt.tight_layout(pad=0.1)
fig.savefig(f"{OUT}/Fig2B_pie.tiff", dpi=200, bbox_inches="tight", pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)

# 범례 단독 (가로 띠, 6.02:1)
fig = plt.figure(figsize=(6.02, 1.0))
fig.legend(handles=[Patch(facecolor=NPG[i], label=f"{LABEL[i]} ({pct[i]:.1f}%)") for i in range(6)],
           loc="center", ncol=3, fontsize=12, frameon=True, handlelength=1.4, columnspacing=1.6)
# bbox_inches를 지정하지 않아 6.02:1 캔버스 비율을 그대로 유지 (pptx 도형 1.83x0.30in에 맞춤)
fig.savefig(f"{OUT}/Fig2B_legend.tiff", dpi=300, pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)

# Fig2A TIFF (pptx 삽입용, 원본과 동일 종횡비 유지)
from PIL import Image as _Im
_Im.open(f"{OUT}/Fig2A_SSR_type.png").convert("RGB").save(f"{OUT}/Fig2A_SSR_type.tiff", compression="tiff_lzw")
for f in ("Fig2A_SSR_type.tiff","Fig2B_pie.tiff","Fig2B_legend.tiff"):
    im=_Im.open(f"{OUT}/{f}"); print(f"  {f}: {im.size} 비율 {im.size[0]/im.size[1]:.3f}")
