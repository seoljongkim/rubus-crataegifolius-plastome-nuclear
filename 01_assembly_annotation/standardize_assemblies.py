#!/usr/bin/env python3
"""
Standardize 10 GetOrganelle chloroplast assemblies against the DAHO reference.

Step 1  – PATH selection (SSC orientation)
    Slide 60 bp anchors through the DAHO SSC region.
    Select the PATH whose SSC matches DAHO orientation.

Step 2  – Circular rotation
    Slide 40 bp windows from ref[0] outward (step 5 bp, up to ref[2000]).
    Use the FIRST window found in the assembly.
    rot_pos = (V - R) % L  →  where ref[R] is at sel_seq[V], so after rotation
    std[0] ≈ ref[0]  (the assembly is anchored to DAHO's start region).

Validation
    Report the signed positional offset (delta = landing − expected position)
    for 10 reference anchors spread across the full genome.
    A correctly rotated sequence shows |delta| ≪ genome_size everywhere,
    with gradual drift reflecting the net indel accumulation between varieties.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import os

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

REF_PATH  = f"{BASE}/Rubus_crataegifolius_DAHO.fasta"
GETORG    = f"{BASE}/getorganelle_output"
OUT_DIR   = f"{BASE}/standardized"

# DAHO structure (0-based): LSC 0–85369 | IR1 85370–111157 | SSC 111158–129875 | IR2 129876–155663
SSC_START  = 111158
SSC_END    = 129876

ANCHOR_LEN = 40    # short enough to avoid spanning an indel
SSC_STEP   = 10    # step for SSC discrimination search
ROT_STEP   = 5     # fine step for rotation anchor (ref[0:2000])
ROT_RANGE  = 2000  # search within first 2000 bp of reference for rotation anchor

VARIETIES = [
    "DAHK_2", "DAH_B1", "DB_A1",  "DS_A1",  "GB048",
    "GB066",  "GDCH_1", "GDM_A1", "GDW_5",  "HGB_A1",
]

os.makedirs(OUT_DIR, exist_ok=True)

# ── helpers ───────────────────────────────────────────────────────────────────
def find_circ(seq: str, pat: str) -> int:
    pos = (seq + seq).find(pat)
    return pos if 0 <= pos < len(seq) else -1

def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def rotate(seq: str, pos: int) -> str:
    return seq[pos:] + seq[:pos]

# ── load reference ────────────────────────────────────────────────────────────
ref_rec = next(SeqIO.parse(REF_PATH, "fasta"))
ref_seq = str(ref_rec.seq).upper()
ref_len = len(ref_seq)

print(f"Reference : {ref_rec.id}  ({ref_len:,} bp)")
print(f"SSC range : {SSC_START+1:,}–{SSC_END:,}  (0-based)")
print(f"Rotation  : first-found {ANCHOR_LEN} bp anchor in ref[0:{ROT_RANGE}], step {ROT_STEP}")
print()

# ── step 1: SSC orientation ───────────────────────────────────────────────────
def select_path_by_ssc(p1: str, p2: str) -> tuple:
    mid  = (SSC_START + SSC_END) // 2
    half = (SSC_END - SSC_START) // 2
    for R in range(mid - half//2, mid + half//2 - 60, SSC_STEP):
        anc = ref_seq[R:R + 60]
        f1  = find_circ(p1, anc)
        f2  = find_circ(p2, anc)
        if f1 != -1 and f2 == -1:
            return p1, "PATH1", R
        if f2 != -1 and f1 == -1:
            return p2, "PATH2", R
    return p1, "PATH1", -1   # fallback

# ── step 2: rotation using first-found anchor near ref[0] ────────────────────
def find_rotation(seq: str) -> tuple:
    """
    Slide ANCHOR_LEN windows from ref[0] outward.
    Return (rot_pos, ref_anchor_start, direction) or raise ValueError.
    """
    L = len(seq)
    for R in range(0, ROT_RANGE - ANCHOR_LEN + 1, ROT_STEP):
        anc = ref_seq[R:R + ANCHOR_LEN]
        V   = find_circ(seq, anc)
        if V != -1:
            return (V - R) % L, R, "fwd"

    # safety: try RC of the whole sequence
    seq_rc = rc(seq)
    for R in range(0, ROT_RANGE - ANCHOR_LEN + 1, ROT_STEP):
        anc = ref_seq[R:R + ANCHOR_LEN]
        V   = find_circ(seq_rc, anc)
        if V != -1:
            return (V - R) % L, R, "rc"

    raise ValueError("No rotation anchor found")

# ── main loop ─────────────────────────────────────────────────────────────────
results = []

for var in VARIETIES:
    base  = f"{GETORG}/{var}/embplant_pt.K115.complete.graph1"
    p1    = str(next(SeqIO.parse(f"{base}.1.path_sequence.fasta", "fasta")).seq).upper()
    p2    = str(next(SeqIO.parse(f"{base}.2.path_sequence.fasta", "fasta")).seq).upper()

    sel_seq, sel_label, ssc_ref_pos = select_path_by_ssc(p1, p2)

    try:
        rot_pos, rot_ref_start, direction = find_rotation(sel_seq)
    except ValueError:
        print(f"  {var}: ERROR — no rotation anchor found! Skipping.")
        continue

    if direction == "rc":
        sel_seq    = rc(sel_seq)
        sel_label += "(RC)"

    rotated = rotate(sel_seq, rot_pos)

    # write
    out_path = f"{OUT_DIR}/{var}_std.fasta"
    with open(out_path, "w") as fh:
        fh.write(f">{var}_std\n")
        for i in range(0, len(rotated), 60):
            fh.write(rotated[i:i+60] + "\n")

    results.append(dict(
        variety       = var,
        selected      = sel_label,
        ssc_ref_pos   = ssc_ref_pos,
        rot_pos       = rot_pos,
        rot_ref_start = rot_ref_start,
        direction     = direction,
        length        = len(rotated),
    ))

# ── validation: delta pattern across genome ───────────────────────────────────
test_positions = [0, 1000, 5000, 15000, 30000, 50000, 80000, 111000, 130000, 150000]

print(f"{'='*90}")
print("Rotation summary and positional delta check")
print(f"  delta = (anchor_landing − expected_position)  [signed, bp]")
print(f"  A correct rotation shows small |delta| with gradual drift due to indel accumulation.")
print(f"{'='*90}")

# header
hdr = f"{'Variety':<10} {'Path':<11} {'RotPos':>8}  {'Anchor@':>7}  {'Len':>8}  "
for p in test_positions:
    hdr += f" {'ref'+str(p//1000)+'k':>7}"
print(hdr)
print("-" * len(hdr))

for r in results:
    std_seq = str(next(SeqIO.parse(f"{OUT_DIR}/{r['variety']}_std.fasta", "fasta")).seq).upper()
    L = len(std_seq)

    row = (f"{r['variety']:<10} {r['selected']:<11} {r['rot_pos']:>8,}  "
           f"ref{r['rot_ref_start']:>4}  {r['length']:>8,}  ")
    for Rp in test_positions:
        anc = ref_seq[Rp:Rp + 40]
        pos = (std_seq + std_seq).find(anc)
        if pos == -1 or pos >= L:
            row += f" {'N/F':>7}"
        else:
            delta = (pos - Rp + L//2) % L - L//2
            row += f" {delta:>+7}"
    print(row)

print(f"\nFiles written : {len(results)}/{len(VARIETIES)}")
print(f"Output        : {OUT_DIR}/")

# ── SSC orientation confirmation ──────────────────────────────────────────────
print(f"\n{'='*60}")
print("SSC orientation and path selection:")
for r in results:
    ssc_tag = "discriminating SSC anchor" if r['ssc_ref_pos'] != -1 else "fallback (PATH1)"
    print(f"  {r['variety']:<10}: {r['selected']:<11}  ({ssc_tag} @ ref {r['ssc_ref_pos']:,})")
