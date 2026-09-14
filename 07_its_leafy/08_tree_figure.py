#!/usr/bin/env python3
"""08_tree_figure.py — ITS ML 계통수 그림 (SH-aLRT/UFBoot 표기)"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
import os

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
os.chdir(f"{BASE}/its_analysis")

t = Phylo.read("ITS_ML2.treefile", "newick")
t.root_with_outgroup("Rosa_multiflora_MW309407")
t.ladderize(reverse=False)

CULT = {"DAHE_Dongakheuk":"'Dongakheuk' (DAHE)","DAHO_Dongakhong":"'Dongakhong' (DAHO)","DAHO_Dongakong":"'Dongakhong' (DAHO)",
        "DB_Dabok":"'Dabok' (DB)","DS_Daesung":"'Daesung' (DS)","GDC_Geumdongcheon":"'Geumdongcheon' (GDC)",
        "GDM_Goldmoon":"'Goldmoon' (GDM)","GDW_Geumdongwang":"'Geumdongwang' (GDW)",
        "HGB_Hwangeumboll":"'Hwangeumboll' (HGB)"}
def pretty(n):
    if n in CULT: return CULT[n]
    p = n.split("_")
    if n.startswith("Rcrataegifolius"): return f"R. crataegifolius {p[1]} ({p[-1]})"
    if n.startswith("Rtakesimensis"):   return f"R. takesimensis {p[1]} ({p[-1]})"
    return f"{p[0][0]}. {p[1]} ({p[-1]})"
for c in t.get_terminals(): c.name = pretty(c.name)

sup = {}
for c in t.get_nonterminals():
    s = c.confidence if isinstance(c.confidence, str) else (c.name or "")
    c.name = None
    if s and "/" in s:
        try:
            sh, bb = (float(x) for x in s.split("/"))
            if bb >= 50 or sh >= 70: sup[c] = f"{sh:.0f}/{bb:.0f}"
        except ValueError: pass

RED = {"'Daesung' (DS)", "'Hwangeumboll' (HGB)"}
fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(t, axes=ax, do_show=False, show_confidence=False,
           branch_labels=lambda c: sup.get(c, ""),
           label_func=lambda c: c.name if c.is_terminal() else None,
           label_colors=lambda n: ("#c0392b" if n in RED
                                   else "#1f6f8b" if n.startswith("R. takesimensis") else "black"))
ax.set_xlabel("substitutions per site", fontsize=10)
ax.set_ylabel(""); ax.set_yticks([])
ax.set_title("nrDNA ITS1–5.8S–ITS2 maximum-likelihood phylogeny (TN+G4, 627 bp)\n"
             "Node labels: SH-aLRT / UFBoot (%)  ·  1,000 replicates each",
             fontsize=10.5)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig("Fig6_ITS_ML_tree.png", dpi=300)
plt.savefig("Fig6_ITS_ML_tree.pdf")
print(f"저장 완료 · 지지도 표시 {len(sup)}개 노드")
