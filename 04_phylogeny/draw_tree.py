#!/usr/bin/env python3
"""draw_tree.py — Fig. 4: 전장 엽록체 기반 ML 계통수 (21 taxa)
GB048/GB066 제외 · 'Dongakhong' 표기 · 외군은 P. mume.
외군 가지가 내군 분기의 20배라 내군이 뭉개지므로, 외군 가지만 축약하고 이중 사선(//)으로 표시."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
import os, re

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

os.chdir(f"{BASE}/phylogeny21")
OUT = "Prunus_mume_MW759300"
t = Phylo.read("ML21.treefile", "newick")
t.root_with_outgroup(OUT)
t.ladderize(reverse=True)

# 루트 자식 중 외군이 아닌 쪽(=내군 줄기)이 긴 가지다
ing = [c for c in t.root.clades if OUT not in [x.name for x in c.get_terminals()]]
ing = ing[0] if ing else t.root.clades[0]
TRUE_LEN = ing.branch_length or 0.0
ing.branch_length = 0.008                      # 내군 줄기 축약

CULT = {"Rc_Dongakheuk":"Dongakheuk","Rc_Dongakhong":"Dongakhong","Rc_Dabok":"Dabok",
        "Rc_Daesung":"Daesung","Rc_Geumdongcheon":"Geumdongcheon","Rc_Goldmoon":"Goldmoon",
        "Rc_Geumdongwang":"Geumdongwang","Rc_Hwangeumboll":"Hwangeumboll"}
IT = lambda s: "$\\it{" + s.replace(" ", "\\ ") + "}$"

def label(n):
    if n in CULT: return IT("R. crataegifolius") + f" ‘{CULT[n]}’"
    m = re.match(r"([A-Z][a-z]+)_([a-z]+)_(.+)", n)
    g, sp, acc = m.group(1), m.group(2), m.group(3).replace("_", "")
    return IT(("P." if g == "Prunus" else "R.") + " " + sp) + f" {acc}"

for c in t.get_terminals(): c.name = label(c.name)
sup = {}
for c in t.get_nonterminals():
    s = c.confidence if isinstance(c.confidence, str) else (c.name or "")
    if not s and c.confidence is not None: s = str(c.confidence)
    c.name = None
    if s and "/" in s:
        sh, bb = s.split("/")
        if float(bb) >= 50: sup[c] = f"{float(sh):.0f}/{float(bb):.0f}"

plt.rcParams.update({"font.family":"serif","font.size":10.5})
fig, ax = plt.subplots(figsize=(11.5, 7.8))
Phylo.draw(t, axes=ax, do_show=False, show_confidence=False,
           branch_labels=lambda c: sup.get(c, ""),
           label_func=lambda c: c.name if c.is_terminal() else None)
# 지지도 라벨: 노드 왼쪽에 붙여 말단 라벨과의 겹침 방지
for txt in ax.texts:
    s = txt.get_text()
    if re.fullmatch(r"\d+/\d+", s.strip()):
        x, y = txt.get_position()
        txt.set_position((x, y)); txt.set_ha("right"); txt.set_va("bottom")
        txt.set_fontsize(7.6); txt.set_color("#333333")
# 내군 줄기(0 → 0.008 구간의 수평선)를 찾아 그 위에 절단 표시
stem_y = None
segs = []
for coll in ax.collections:
    try: segs += list(coll.get_segments())
    except Exception: pass
for sg in segs:
    (x0, y0), (x1, y1) = sg[0], sg[-1]
    if abs(y0-y1) < 1e-9 and abs((x1-x0)-0.008) < 5e-4 and x0 < 1e-3:
        stem_y = y0; break
if stem_y is not None:
    for x0 in (0.0030, 0.0045):
        ax.plot([x0, x0+0.0013], [stem_y-0.35, stem_y+0.35], color="k", lw=1.3, zorder=5)
    ax.text(0.0038, stem_y+0.55, "//", ha="center", va="bottom", fontsize=1, color="white")
else:
    print("  [경고] 내군 줄기를 찾지 못해 절단 표시 생략")
ax.set_xlabel("substitutions per site", fontsize=11.5)
ax.set_ylabel(""); ax.set_yticks([])
for s in ("top","right","left"): ax.spines[s].set_visible(False)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(f"{BASE}/figures_revision/Fig4_cp_phylogeny.{ext}",
                dpi=300 if ext=="png" else None, bbox_inches="tight")
print(f"저장 완료 · 말단 {len(t.get_terminals())}, 지지도 {len(sup)} | 외군 실제 가지길이 {TRUE_LEN:.4f} → 0.010으로 축약")
