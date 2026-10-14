#!/usr/bin/env python3
"""generate_nuclear_figure.py — Fig. 6: 핵 SNP 기반 품종 관계 (HEB 리뷰어 R1-b, R2-2 대응)
(A) 핵 1-IBS 거리행렬  (B) 엽록체 거리 vs 핵 거리  (C) 핵 거리 UPGMA 덴드로그램"""
import numpy as np, itertools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
import numpy as _np

import os

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

K=f"{BASE}/kaks_output"; OUT=f"{BASE}/figures_revision"
VARS=["DAHE","DAHO","DB","DS","GDC","GDM","GDW","HGB"]
FULL={"DAHE":"Dongakheuk","DAHO":"Dongakhong","DB":"Dabok","DS":"Daesung",
      "GDC":"Geumdongcheon","GDM":"Goldmoon","GDW":"Geumdongwang","HGB":"Hwangeumboll"}
GROUP={"DAHE":"A","DAHO":"A","GDW":"A","DB":"B","GDC":"B","DS":"–","GDM":"–","HGB":"–"}
NPG=["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4"]

lab=open(f"{K}/nuc_dist_labels.txt").read().split()
nuc=np.load(f"{K}/nuc_dist_poly.npy")   # 품종 간 다형 4,872,213개 기준 (종간 고정차이 제외)
cp=np.load(f"{K}/cp_dist.npy")

plt.rcParams.update({"font.family":"serif","font.size":11})
fig=plt.figure(figsize=(17.5,6.4))
gs=fig.add_gridspec(1,3,width_ratios=[1.05,1.0,1.05],wspace=0.42,
                    left=.045,right=.985,top=.90,bottom=.24)

# ---------- (A) 거리행렬 히트맵 ----------
ax=fig.add_subplot(gs[0])
cmap=LinearSegmentedColormap.from_list("nb",["#FFFFFF","#9EC7E0","#3C5488"])
m=np.ma.masked_where(np.eye(8,dtype=bool),nuc)
im=ax.imshow(m,cmap=cmap,vmin=nuc[nuc>0].min(),vmax=nuc.max())
for i in range(8):
    for j in range(8):
        if i==j: continue
        v=nuc[i,j]
        ax.text(j,i,f"{v:.3f}",ha="center",va="center",fontsize=10.5,
                color="white" if v>nuc.max()*0.82 else "black")
    ax.text(i,i,"—",ha="center",va="center",color="#999999")
ax.set_xticks(range(8)); ax.set_yticks(range(8))
ax.set_xticklabels([FULL[v] for v in VARS],rotation=45,ha="right",fontsize=13)
ax.set_yticklabels([FULL[v] for v in VARS],fontsize=13)
# cp 동일 그룹 쌍에 테두리
for i,j in itertools.combinations(range(8),2):
    if cp[i,j]==0:
        for (a,b) in ((i,j),(j,i)):
            ax.add_patch(plt.Rectangle((b-0.5,a-0.5),1,1,fill=False,edgecolor="#E64B35",lw=2.0))
fig.colorbar(im,ax=ax,fraction=0.045,pad=0.03,label="nuclear 1–IBS distance")
from matplotlib.patches import Rectangle as _R
ax.legend(handles=[_R((0,0),1,1,fill=False,edgecolor="#E64B35",lw=2.0,
                      label="identical cp genomes (0 SNPs)")],
          fontsize=12.5,loc="upper center",bbox_to_anchor=(.5,-.30),frameon=False,handletextpad=.6)
ax.set_title("(a)  Nuclear pairwise distance",loc="left",fontsize=14.5,fontweight="bold",pad=12)

# ---------- (B) cp vs 핵 산점도 ----------
ax=fig.add_subplot(gs[1])
iu=np.triu_indices(8,1); x=cp[iu]; y=nuc[iu]
same=x==0
ax.scatter(x[~same],y[~same],s=62,color=NPG[3],alpha=.82,edgecolor="white",lw=.8,
           label="different cp genomes",zorder=3)
ax.scatter(x[same],y[same],s=110,color=NPG[0],edgecolor="black",lw=1.0,marker="D",
           label="identical cp genomes (0 SNPs)",zorder=4)
# cp 동일 쌍 주석: 풀네임이라 라벨이 길어져, 거리순으로 위·아래 번갈아 배치해 겹침을 피한다
_id_pairs=sorted([(a,b) for a,b in zip(*iu) if cp[a,b]==0],key=lambda t:nuc[t])
for _k,(a,b) in enumerate(_id_pairs):
    ax.annotate(f"{FULL[VARS[a]]}–{FULL[VARS[b]]}",(cp[a,b],nuc[a,b]),
                textcoords="offset points",xytext=(11,7 if _k%2 else -12),
                fontsize=10.5,color=NPG[0])
ax.set_ylim(ax.get_ylim()[0]-0.008,ax.get_ylim()[1])   # 아래쪽 라벨이 축선에 닿지 않도록
rho,_=spearmanr(x,y)
_rng=_np.random.default_rng(20260818); _N=99999; _c=0
for _ in range(_N):
    _pm=_rng.permutation(8)
    if abs(spearmanr(cp[_np.ix_(_pm,_pm)][iu], y)[0])>=abs(rho): _c+=1
p=(_c+1)/(_N+1)   # Mantel 순열검정 (쌍별 거리의 비독립성 보정)
ax.set_xlabel("chloroplast pairwise SNPs",fontsize=13)
ax.set_ylabel("nuclear 1–IBS distance",fontsize=13)
ax.set_title("(b)  Chloroplast vs. nuclear divergence",loc="left",fontsize=14.5,fontweight="bold",pad=12)
ax.text(.975,.50,f"Mantel ρ = {rho:.2f}\nP = {p:.3f}\n99,999 permutations",transform=ax.transAxes,
        ha="right",va="center",fontsize=12,
        bbox=dict(boxstyle="round,pad=0.42",fc="white",ec="#BBBBBB"))
ax.legend(fontsize=12.5,loc="upper center",bbox_to_anchor=(.5,-.155),ncol=1,frameon=False,
          handletextpad=.6,labelspacing=.5)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.grid(alpha=.22,linestyle=":")

# ---------- (C) UPGMA 덴드로그램 ----------
ax=fig.add_subplot(gs[2])
Z=linkage(squareform(nuc,checks=False),method="average")
GC={"A":NPG[0],"B":NPG[2],"–":"#555555"}
d=dendrogram(Z,labels=[f"{FULL[v]}" for v in VARS],orientation="left",ax=ax,
             color_threshold=0,above_threshold_color="#333333")
for t in ax.get_yticklabels():
    v=[k for k in VARS if FULL[k]==t.get_text()][0]
    t.set_color(GC[GROUP[v]]); t.set_fontsize(12.5)
ax.set_xlabel("nuclear 1–IBS distance",fontsize=13)
ax.set_title("(c)  UPGMA clustering (nuclear SNPs)",loc="left",fontsize=14.5,fontweight="bold",pad=12)
ax.yaxis.tick_right()
for s in ("top","right","left"): ax.spines[s].set_visible(False)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],color=NPG[0],lw=3,label="cp haplotype group A"),
                   Line2D([],[],color=NPG[2],lw=3,label="cp haplotype group B"),
                   Line2D([],[],color="#555555",lw=3,label="unique cp haplotype")],
          fontsize=12.5,loc="upper center",bbox_to_anchor=(.5,-.155),ncol=1,frameon=False,
          handletextpad=.6,labelspacing=.5)

fig.savefig(f"{OUT}/Fig6_nuclear_SNP.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{OUT}/Fig6_nuclear_SNP.pdf",bbox_inches="tight")
print(f"저장: {OUT}/Fig6_nuclear_SNP.png / .pdf")
print(f"Mantel ρ={rho:.3f}, P={p:.4f}")
print(f"핵 거리 범위 {y.min():.5f}~{y.max():.5f}, cp 동일 4쌍 평균 {y[same].mean():.5f}")
