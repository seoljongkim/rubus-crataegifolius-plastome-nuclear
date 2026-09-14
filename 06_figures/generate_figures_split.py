#!/usr/bin/env python3
"""Split panel figures into individual files (8 varieties, GB048/GB066 excluded)."""

import os, csv, re, sys
csv.field_size_limit(sys.maxsize)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

OUT_DIR = f'{BASE}/figures'
os.makedirs(OUT_DIR, exist_ok=True)

VARIETIES = ['DAHE', 'DAHO', 'DB', 'DS', 'GDC', 'GDM', 'GDW', 'HGB']
ALL_10 = ['DAHE', 'DAHO', 'DB', 'DS', 'GB048', 'GB066', 'GDC', 'GDM', 'GDW', 'HGB']

COLORS = {
    'mono': '#E53935', 'di': '#FB8C00', 'tri': '#43A047', 'tetra': '#1E88E5',
    'penta': '#8E24AA', 'hexa': '#6D4C41',
    'Forward': '#1976D2', 'Palindromic': '#D32F2F', 'Reverse': '#388E3C', 'Complement': '#F57C00',
    'SNP': '#1565C0', 'InDel': '#C62828', 'Complex': '#F9A825'
}

type_labels = {1:'Mono', 2:'Di', 3:'Tri', 4:'Tetra', 5:'Penta', 6:'Hexa'}
type_colors = [COLORS['mono'], COLORS['di'], COLORS['tri'],
               COLORS['tetra'], COLORS['penta'], COLORS['hexa']]


def save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, f'{name}.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  -> {name}.png/pdf")


def parse_ssr():
    counts = {}
    for v in VARIETIES:
        counts[v] = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0}
        with open(f'{BASE}/misa_output/{v}_std.fasta.statistics') as f:
            in_dist = False
            for line in f:
                line = line.strip()
                if line.startswith('Unit size\tNumber'):
                    in_dist = True; continue
                if in_dist and line:
                    parts = line.split('\t')
                    if len(parts) == 2 and parts[0].isdigit():
                        counts[v][int(parts[0])] = int(parts[1])
                    elif not line[0].isdigit():
                        in_dist = False
    return counts


# ===================== Fig 3A: SSR Type Distribution =====================
def fig3a():
    print("Fig 3A: SSR Type Distribution...")
    counts = parse_ssr()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(VARIETIES))
    bottom = np.zeros(len(VARIETIES))
    for i, us in enumerate([1,2,3,4,5,6]):
        vals = [counts[v][us] for v in VARIETIES]
        ax.bar(x, vals, 0.6, bottom=bottom, label=type_labels[us],
               color=type_colors[i], edgecolor='white', linewidth=0.5)
        bottom += np.array(vals)
    for i in range(len(VARIETIES)):
        ax.text(i, bottom[i]+5, str(int(bottom[i])), ha='center', va='bottom',
                fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of SSRs', fontsize=12)
    ax.set_title('SSR Type Distribution', fontsize=14)
    ax.legend(title='Repeat Unit', fontsize=10, title_fontsize=11, loc='upper right')
    plt.tight_layout()
    save(fig, 'Fig3A_SSR_type')


# ===================== Fig 3B: SSR Type Proportion =====================
def fig3b():
    print("Fig 3B: SSR Type Proportion...")
    counts = parse_ssr()
    fig, ax = plt.subplots(figsize=(8, 7))
    avg = [np.mean([counts[v][us] for v in VARIETIES]) for us in [1,2,3,4,5,6]]
    labels_pie = [f'{type_labels[i+1]}\n({int(s)})' for i, s in enumerate(avg)]
    wedges, texts, autotexts = ax.pie(avg, labels=labels_pie, colors=type_colors,
                                       autopct='%1.1f%%', startangle=90,
                                       pctdistance=0.75, textprops={'fontsize': 10})
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title('Average SSR Type Proportion', fontsize=14)
    plt.tight_layout()
    save(fig, 'Fig3B_SSR_proportion')


# ===================== Fig 4A: Tandem Repeat Length =====================
def fig4a():
    print("Fig 4A: Tandem Repeat Length Distribution...")
    all_len = []
    for v in VARIETIES:
        with open(f'{BASE}/trf_output/{v}_std.fasta.2.7.7.80.10.50.500.dat', 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
        for line in content.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0].isdigit():
                try:
                    all_len.append(int(parts[1]) - int(parts[0]) + 1)
                except ValueError:
                    pass

    fig, ax = plt.subplots(figsize=(8, 6))
    bins = [0, 20, 40, 60, 80, 100, 150, 200, 500]
    bin_labels = ['<20', '20-39', '40-59', '60-79', '80-99', '100-149', '150-199', '≥200']
    cnts, _ = np.histogram(all_len, bins=bins)
    ax.bar(range(len(bin_labels)), cnts, color='#5C6BC0', edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(len(bin_labels)))
    ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Repeats', fontsize=12)
    ax.set_xlabel('Repeat Length (bp)', fontsize=12)
    ax.set_title('Tandem Repeat Length Distribution', fontsize=14)
    for i, c in enumerate(cnts):
        if c > 0:
            ax.text(i, c+1, str(c), ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    save(fig, 'Fig4A_Tandem_length')


# ===================== Fig 4B: Dispersed Repeat Type =====================
def parse_dispersed():
    disp_type = {}; disp_size = {}
    with open(f'{BASE}/reputer_output/dispersed_repeat_summary.tsv') as f:
        content = f.read()
    sections = content.split('===')
    for i, sec in enumerate(sections):
        if 'Counts by Type' in sec:
            for line in sections[i+1].strip().split('\n'):
                parts = line.strip().split('\t')
                if len(parts) >= 7 and parts[0] in VARIETIES:
                    disp_type[parts[0]] = {'Forward':int(parts[3]), 'Palindromic':int(parts[4]),
                                           'Reverse':int(parts[5]), 'Complement':int(parts[6])}
        elif 'Size Distribution of Dispersed' in sec:
            for line in sections[i+1].strip().split('\n'):
                parts = line.strip().split('\t')
                if len(parts) >= 5 and parts[0] in VARIETIES:
                    disp_size[parts[0]] = {'30-44':int(parts[1]), '45-59':int(parts[2]),
                                           '60-89':int(parts[3]), '≥90':int(parts[4])}
    return disp_type, disp_size


def fig4b():
    print("Fig 4B: Dispersed Repeat Type Distribution...")
    disp_type, _ = parse_dispersed()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(VARIETIES)); w = 0.2
    for i, t in enumerate(['Forward', 'Palindromic', 'Reverse', 'Complement']):
        vals = [disp_type.get(v,{}).get(t,0) for v in VARIETIES]
        if max(vals) > 0:
            col = [COLORS['Forward'], COLORS['Palindromic'], COLORS['Reverse'], COLORS['Complement']][i]
            ax.bar(x + i*w - 1.5*w, vals, w, label=t, color=col, edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Repeats', fontsize=12)
    ax.set_title('Dispersed Repeat Type Distribution', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    plt.tight_layout()
    save(fig, 'Fig4B_Dispersed_type')


# ===================== Fig 4C: Dispersed Repeat Size =====================
def fig4c():
    print("Fig 4C: Dispersed Repeat Size Distribution...")
    _, disp_size = parse_dispersed()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(VARIETIES))
    sr_list = ['30-44', '45-59', '60-89', '≥90']
    sr_cols = ['#81D4FA', '#4FC3F7', '#29B6F6', '#0288D1']
    bottom = np.zeros(len(VARIETIES))
    for i, sr in enumerate(sr_list):
        vals = np.array([disp_size.get(v,{}).get(sr,0) for v in VARIETIES])
        ax.bar(x, vals, 0.6, bottom=bottom, label=f'{sr} bp', color=sr_cols[i],
               edgecolor='white', linewidth=0.5)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Repeats', fontsize=12)
    ax.set_title('Dispersed Repeat Size Distribution', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    plt.tight_layout()
    save(fig, 'Fig4C_Dispersed_size')


# ===================== Fig 6 Variant helpers =====================
def parse_variants_8():
    excluded = {'GB048', 'GB066'}
    region_counts = defaultdict(lambda: {'SNP':0, 'InDel':0, 'Complex':0})
    valid_variants = []

    with open(f'{BASE}/variant_output/all_variants.csv') as f:
        for row in csv.DictReader(f):
            details = row['Details']; vtype = row['Type']; region = row['Region']
            if vtype == 'SNP':
                allele_groups = {}
                for part in details.split('; '):
                    if ':' in part:
                        allele, vars_str = part.split(':', 1)
                        vars_8 = [v.strip() for v in vars_str.split(',') if v.strip() not in excluded]
                        if vars_8:
                            allele_groups[allele] = vars_8
                if len(allele_groups) >= 2:
                    valid_variants.append(row)
                    region_counts[region]['SNP'] += 1
            elif vtype == 'InDel':
                m = re.match(r'gap_in=([^;]+)', details)
                if m:
                    gap_8 = [v.strip() for v in m.group(1).split(',') if v.strip() not in excluded]
                    non_gap_8 = [v for v in VARIETIES if v not in gap_8]
                    if gap_8 and non_gap_8:
                        valid_variants.append(row)
                        region_counts[region]['InDel'] += 1
                else:
                    valid_variants.append(row)
                    region_counts[region]['InDel'] += 1
            else:
                valid_variants.append(row)
                region_counts[region]['Complex'] += 1

    # Private variants
    private = {v: {'SNPs':0, 'InDels':0} for v in VARIETIES}
    for row in valid_variants:
        details = row['Details']; vtype = row['Type']
        if vtype == 'SNP':
            for part in details.split('; '):
                if ':' in part:
                    _, vars_str = part.split(':', 1)
                    vl = [v.strip() for v in vars_str.split(',') if v.strip() not in excluded]
                    if len(vl) == 1 and vl[0] in VARIETIES:
                        private[vl[0]]['SNPs'] += 1
        elif vtype == 'InDel':
            m = re.match(r'gap_in=([^;]+)', details)
            if m:
                gap = [v.strip() for v in m.group(1).split(',') if v.strip() not in excluded]
                non_gap = [v for v in VARIETIES if v not in gap]
                if len(gap) == 1 and gap[0] in VARIETIES:
                    private[gap[0]]['InDels'] += 1
                elif len(non_gap) == 1 and non_gap[0] in VARIETIES:
                    private[non_gap[0]]['InDels'] += 1

    # Pairwise SNP from original 10x10
    snp_matrix_10 = {}; mv10 = []
    with open(f'{BASE}/variant_output/variant_summary.txt') as f:
        content = f.read()
    in_m = False
    for line in content.split('\n'):
        if 'Pairwise SNP counts' in line:
            in_m = True; continue
        if in_m:
            parts = line.strip().split()
            if parts and parts[0] in ALL_10:
                vals = []
                for p in parts[1:]:
                    if p == '—': vals.append(0)
                    else:
                        try: vals.append(int(p))
                        except ValueError: pass
                snp_matrix_10[parts[0]] = vals
                mv10.append(parts[0])

    n8 = len(VARIETIES)
    mat = np.zeros((n8, n8))
    for i, vi in enumerate(VARIETIES):
        if vi in snp_matrix_10:
            rv = snp_matrix_10[vi]
            for j, vj in enumerate(VARIETIES):
                j10 = mv10.index(vj) if vj in mv10 else -1
                if 0 <= j10 < len(rv):
                    mat[i][j] = rv[j10]

    return region_counts, private, mat


# ===================== Fig 6A: Variants by Region =====================
def fig6a():
    print("Fig 6A: Variants by Genome Region...")
    region_counts, _, _ = parse_variants_8()
    region_sizes = {'LSC': 85374, 'IRb': 25786, 'SSC': 18718, 'IRa': 25786}
    regions_order = ['LSC', 'IRb', 'SSC', 'IRa']

    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(4); w = 0.25
    snps = [region_counts[r]['SNP'] for r in regions_order]
    indels = [region_counts[r]['InDel'] for r in regions_order]
    cmplx = [region_counts[r]['Complex'] for r in regions_order]

    ax.bar(x-w, snps, w, label='SNPs', color=COLORS['SNP'])
    ax.bar(x, indels, w, label='InDels', color=COLORS['InDel'])
    ax.bar(x+w, cmplx, w, label='Complex', color=COLORS['Complex'])
    ax.set_xticks(x)
    ax.set_xticklabels(regions_order, fontsize=11)
    ax.set_ylabel('Number of Variants', fontsize=12)
    ax.set_title('Variants by Genome Region', fontsize=14)
    ax.legend(fontsize=10)

    for i, r in enumerate(regions_order):
        total = region_counts[r]['SNP'] + region_counts[r]['InDel'] + region_counts[r]['Complex']
        per_kb = total / (region_sizes[r] / 1000)
        ax.text(i, total+5, f'{per_kb:.2f}/kb', ha='center', va='bottom', fontsize=9, color='gray')

    plt.tight_layout()
    save(fig, 'Fig6A_Variant_region')


# ===================== Fig 6B: Private Variants =====================
def fig6b():
    print("Fig 6B: Private Variants per Variety...")
    _, private, _ = parse_variants_8()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(len(VARIETIES)); w = 0.35
    ps = [private[v]['SNPs'] for v in VARIETIES]
    pi = [private[v]['InDels'] for v in VARIETIES]
    ax.bar(x - w/2, ps, w, label='Private SNPs', color=COLORS['SNP'])
    ax.bar(x + w/2, pi, w, label='Private InDels', color=COLORS['InDel'])
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Private Variants', fontsize=12)
    ax.set_title('Private Variants per Variety', fontsize=14)
    ax.legend(fontsize=10)
    plt.tight_layout()
    save(fig, 'Fig6B_Variant_private')


# ===================== Fig 6C: Pairwise SNP Heatmap =====================
def fig6c():
    print("Fig 6C: Pairwise SNP Heatmap...")
    _, _, mat = parse_variants_8()
    n8 = len(VARIETIES)
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, cmap='YlOrRd', aspect='equal')
    ax.set_xticks(range(n8))
    ax.set_yticks(range(n8))
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(VARIETIES, fontsize=10)
    ax.set_title('Pairwise SNP Counts', fontsize=14)
    for i in range(n8):
        for j in range(n8):
            val = int(mat[i][j])
            color = 'white' if val > 200 else 'black'
            ax.text(j, i, str(val), ha='center', va='center', fontsize=8, color=color)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('SNP Count', fontsize=11)
    plt.tight_layout()
    save(fig, 'Fig6C_Variant_pairwise')


if __name__ == '__main__':
    print("=" * 60)
    print("Generating split panel figures (8 varieties)")
    print("=" * 60)
    fig3a(); fig3b()
    fig4a(); fig4b(); fig4c()
    fig6a(); fig6b(); fig6c()
    print("\n" + "=" * 60)
    print(f"All saved to: {OUT_DIR}")
    print("=" * 60)
