#!/usr/bin/env python3
"""
Generate publication figures for Rubus crataegifolius chloroplast genome paper.
8 varieties only (GB048, GB066 excluded).
Fig 2: RSCU bar chart
Fig 3: SSR analysis (type distribution + motif proportion)
Fig 4: Repeat analysis (tandem + dispersed)
Fig 5: IR boundary comparison (IRscope-style)
Fig 6: Structural variant analysis (region + private + pairwise heatmap)
"""

import os, csv, re, sys
csv.field_size_limit(sys.maxsize)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from collections import defaultdict

OUT_DIR = '/home/spenkim/rubus_cp/figures'
os.makedirs(OUT_DIR, exist_ok=True)

VARIETIES = ['DAHE', 'DAHO', 'DB', 'DS', 'GDC', 'GDM', 'GDW', 'HGB']
ALL_10 = ['DAHE', 'DAHO', 'DB', 'DS', 'GB048', 'GB066', 'GDC', 'GDM', 'GDW', 'HGB']
VARIETY_LABELS = {
    'DAHE': 'DAHK_2', 'DAHO': 'DAH_B1', 'DB': 'DB_A1', 'DS': 'DS_A1',
    'GDC': 'GDCH_1', 'GDM': 'GDM_A1', 'GDW': 'GDW_5', 'HGB': 'HGB_A1'
}

# 품종명 (그림 표시용)
CULTIVAR = {
    'DAHE': 'Dongakheuk', 'DAHO': 'Dongakhong', 'DB': 'Dabok', 'DS': 'Daesung',
    'GDC': 'Geumdongcheon', 'GDM': 'Goldmoon', 'GDW': 'Geumdongwang', 'HGB': 'Hwangeumboll'
}

COLORS = {
    'mono': '#E53935', 'di': '#FB8C00', 'tri': '#43A047', 'tetra': '#1E88E5',
    'penta': '#8E24AA', 'hexa': '#6D4C41',
    'Forward': '#1976D2', 'Palindromic': '#D32F2F', 'Reverse': '#388E3C', 'Complement': '#F57C00',
    'SNP': '#1565C0', 'InDel': '#C62828', 'Complex': '#F9A825'
}


###############################################################################
# Fig 2: RSCU bar chart (DAHO representative)
###############################################################################
def fig2_rscu():
    print("Generating Fig 2: RSCU bar chart...")
    codons, amino_acids, data = [], [], []
    with open('/home/spenkim/rubus_cp/rscu_output/rscu_all_varieties.csv') as f:
        for row in csv.DictReader(f):
            if row['AminoAcid'] in ('*', 'Stop'):
                continue
            codons.append(row['Codon'])
            amino_acids.append(row['AminoAcid'])
            data.append(float(row['DAHO']))

    n = len(codons)
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.bar(range(n), data, width=0.8,
           color=['#E53935' if d > 1 else '#1976D2' for d in data],
           edgecolor='white', linewidth=0.3)
    ax.set_xticks(range(n))
    ax.set_xticklabels([f"{c} ({a})" for a, c in zip(amino_acids, codons)],
                       fontsize=8, rotation=60, ha='right', rotation_mode='anchor')
    ax.set_ylabel('RSCU Value', fontsize=12)
    ax.set_title(f"Relative Synonymous Codon Usage (RSCU) - '{CULTIVAR['DAHO']}' (representative)", fontsize=14)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlim(-0.5, n - 0.5)
    ax.legend(handles=[
        mpatches.Patch(color='#E53935', label='RSCU > 1 (Preferred)'),
        mpatches.Patch(color='#1976D2', label='RSCU ≤ 1')
    ], loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'Fig2_RSCU.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'Fig2_RSCU.pdf'), bbox_inches='tight')
    plt.close()
    print("  -> Fig2_RSCU.png/pdf")


###############################################################################
# Fig 3: SSR Analysis
###############################################################################
def fig3_ssr():
    print("Generating Fig 3: SSR analysis...")
    ssr_type_counts = {}
    for v in VARIETIES:
        ssr_type_counts[v] = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0}
        with open(f'/home/spenkim/rubus_cp/misa_output/{v}_std.fasta.statistics') as f:
            in_dist = False
            for line in f:
                line = line.strip()
                if line.startswith('Unit size\tNumber'):
                    in_dist = True; continue
                if in_dist and line:
                    parts = line.split('\t')
                    if len(parts) == 2 and parts[0].isdigit():
                        ssr_type_counts[v][int(parts[0])] = int(parts[1])
                    elif not line[0].isdigit():
                        in_dist = False

    type_labels = {1:'Mono', 2:'Di', 3:'Tri', 4:'Tetra', 5:'Penta', 6:'Hexa'}
    type_colors = [COLORS['mono'], COLORS['di'], COLORS['tri'],
                   COLORS['tetra'], COLORS['penta'], COLORS['hexa']]

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # (A) Stacked bar
    ax = axes[0]
    x = np.arange(len(VARIETIES))
    bottom = np.zeros(len(VARIETIES))
    for i, us in enumerate([1,2,3,4,5,6]):
        vals = [ssr_type_counts[v][us] for v in VARIETIES]
        ax.bar(x, vals, 0.6, bottom=bottom, label=type_labels[us],
               color=type_colors[i], edgecolor='white', linewidth=0.5)
        bottom += np.array(vals)
    for i, v in enumerate(VARIETIES):
        ax.text(i, bottom[i]+5, str(int(bottom[i])), ha='center', va='bottom',
                fontsize=8, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of SSRs', fontsize=11)
    ax.set_title('(A) SSR Type Distribution', fontsize=13)
    ax.legend(title='Repeat Unit', fontsize=9, title_fontsize=10, loc='upper right')

    # (B) Pie chart (average)
    ax = axes[1]
    avg = [np.mean([ssr_type_counts[v][us] for v in VARIETIES]) for us in [1,2,3,4,5,6]]
    labels_pie = [f'{type_labels[i+1]}\n({int(s)})' for i, s in enumerate(avg)]
    wedges, texts, autotexts = ax.pie(avg, labels=labels_pie, colors=type_colors,
                                       autopct='%1.1f%%', startangle=90,
                                       pctdistance=0.75, textprops={'fontsize': 9})
    for t in autotexts:
        t.set_fontsize(8)
    ax.set_title('(B) Average SSR Type Proportion', fontsize=13)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'Fig3_SSR.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'Fig3_SSR.pdf'), bbox_inches='tight')
    plt.close()
    print("  -> Fig3_SSR.png/pdf")


###############################################################################
# Fig 4: Repeat Analysis (Tandem + Dispersed)
###############################################################################
def fig4_repeat():
    print("Generating Fig 4: Repeat analysis...")

    # TRF
    trf_lengths = {}
    for v in VARIETIES:
        trf_lengths[v] = []
        with open(f'/home/spenkim/rubus_cp/trf_output/{v}_std.fasta.2.7.7.80.10.50.500.dat', 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
        for line in content.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0].isdigit():
                try:
                    trf_lengths[v].append(int(parts[1]) - int(parts[0]) + 1)
                except ValueError:
                    pass

    # Dispersed repeats
    disp_type = {}; disp_size = {}
    with open('/home/spenkim/rubus_cp/reputer_output/dispersed_repeat_summary.tsv') as f:
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

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # (A) Tandem length distribution
    ax = axes[0]
    all_len = []
    for v in VARIETIES:
        all_len.extend(trf_lengths[v])
    bins = [0, 20, 40, 60, 80, 100, 150, 200, 500]
    bin_labels = ['<20', '20-39', '40-59', '60-79', '80-99', '100-149', '150-199', '≥200']
    counts, _ = np.histogram(all_len, bins=bins)
    ax.bar(range(len(bin_labels)), counts, color='#5C6BC0', edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(len(bin_labels)))
    ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of Repeats', fontsize=11)
    ax.set_xlabel('Repeat Length (bp)', fontsize=11)
    ax.set_title('(A) Tandem Repeat Length Distribution', fontsize=12)
    for i, c in enumerate(counts):
        if c > 0:
            ax.text(i, c+0.5, str(c), ha='center', va='bottom', fontsize=8)

    # (B) Dispersed type - grouped bar
    ax = axes[1]
    x = np.arange(len(VARIETIES)); w = 0.2
    for i, t in enumerate(['Forward', 'Palindromic', 'Reverse', 'Complement']):
        vals = [disp_type.get(v,{}).get(t,0) for v in VARIETIES]
        if max(vals) > 0:
            col = [COLORS['Forward'], COLORS['Palindromic'], COLORS['Reverse'], COLORS['Complement']][i]
            ax.bar(x + i*w - 1.5*w, vals, w, label=t, color=col, edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of Repeats', fontsize=11)
    ax.set_title('(B) Dispersed Repeat Type Distribution', fontsize=12)
    ax.legend(fontsize=9, loc='upper right')

    # (C) Dispersed size - stacked bar
    ax = axes[2]
    sr_list = ['30-44', '45-59', '60-89', '≥90']
    sr_cols = ['#81D4FA', '#4FC3F7', '#29B6F6', '#0288D1']
    bottom = np.zeros(len(VARIETIES))
    for i, sr in enumerate(sr_list):
        vals = np.array([disp_size.get(v,{}).get(sr,0) for v in VARIETIES])
        ax.bar(x, vals, 0.6, bottom=bottom, label=f'{sr} bp', color=sr_cols[i],
               edgecolor='white', linewidth=0.5)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of Repeats', fontsize=11)
    ax.set_title('(C) Dispersed Repeat Size Distribution', fontsize=12)
    ax.legend(fontsize=9, loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'Fig4_Repeat.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'Fig4_Repeat.pdf'), bbox_inches='tight')
    plt.close()
    print("  -> Fig4_Repeat.png/pdf")


###############################################################################
# Fig 5: IR Boundary Comparison
###############################################################################
def fig5_irscope():
    print("Generating Fig 5: IR boundary comparison...")
    ir_data = {}
    with open('/home/spenkim/rubus_cp/irscope_output/ir_boundary_genes.csv') as f:
        for row in csv.DictReader(f):
            v = row['Variety']
            if v not in VARIETIES:
                continue
            junc = row['Junction']
            ir_data.setdefault(v, {}).setdefault(junc, []).append({
                'gene': row['Gene'], 'pseudo': row.get('Pseudo', ''),
                'dist': int(row['Dist_to_junction']) if row['Dist_to_junction'] else 0,
                'position': row['Position'],
                'junction_pos': int(row['Junction_pos']),
                'length': int(row['Length']), 'strand': row['Strand']
            })

    n_var = len(VARIETIES)
    fig, ax = plt.subplots(figsize=(18, n_var * 1.2 + 4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, n_var * 8 + 10)
    ax.axis('off')

    junctions = ['JLB', 'JSB', 'JSA', 'JLA']
    junction_x = [20, 40, 60, 80]
    junction_names = ['JLB\n(LSC/IRb)', 'JSB\n(IRb/SSC)', 'JSA\n(SSC/IRa)', 'JLA\n(IRa/LSC)']
    region_colors = {'LSC': '#90CAF9', 'IRb': '#FFCC80', 'SSC': '#A5D6A7', 'IRa': '#CE93D8'}

    ax.text(50, n_var*8+7, f'IR Boundary Comparison of {n_var} R. crataegifolius Varieties',
            ha='center', fontsize=16, fontweight='bold')

    header_y = n_var * 8 + 3
    regions = [('LSC',0,20), ('IRb',20,40), ('SSC',40,60), ('IRa',60,80), ('LSC',80,100)]
    for name, x1, x2 in regions:
        ax.add_patch(plt.Rectangle((x1, header_y-1), x2-x1, 2,
                     facecolor=region_colors[name], edgecolor='black', linewidth=1))
        ax.text((x1+x2)/2, header_y, name, ha='center', va='center', fontsize=11, fontweight='bold')

    for jx, jn in zip(junction_x, junction_names):
        ax.axvline(x=jx, ymin=0, ymax=0.9, color='red', linewidth=1.5, linestyle='--', alpha=0.5)
        ax.text(jx, header_y-2.5, jn, ha='center', va='top', fontsize=9, color='red', fontweight='bold')

    for vi, v in enumerate(VARIETIES):
        y_base = (n_var - 1 - vi) * 8 + 2
        label = VARIETY_LABELS.get(v, v)
        ax.text(-1, y_base+1, f"'{CULTIVAR.get(v, v)}'", ha='right', va='center', fontsize=8, fontweight='bold')

        for name, x1, x2 in regions:
            ax.add_patch(plt.Rectangle((x1, y_base), x2-x1, 2,
                         facecolor=region_colors[name], edgecolor='gray', linewidth=0.5, alpha=0.3))

        if v in ir_data:
            for ji, junc in enumerate(junctions):
                jx = junction_x[ji]
                if junc not in ir_data[v]:
                    continue
                genes = [g for g in ir_data[v][junc] if g['gene'] != 'rps12']
                for gi, g in enumerate(genes):
                    gname = f"ψ{g['gene']}" if g['pseudo'] else g['gene']
                    pos = g['position']
                    box_w, box_h = 8, 1.5

                    if pos == 'LSC':
                        bx = jx - box_w - 0.5
                    elif pos == 'spanning':
                        bx = jx - box_w / 2
                    elif pos == 'SSC':
                        bx = jx + 0.5 if junc == 'JSB' else jx - box_w - 0.5
                    elif pos == 'IRb':
                        bx = jx + 0.5
                    elif pos == 'IRa':
                        bx = jx - box_w - 0.5 if junc == 'JLA' else jx + 0.5
                    else:
                        bx = jx + 0.5

                    by = y_base + 2 + 0.3 + gi * 2
                    gc = '#FFEB3B' if g['pseudo'] else '#E3F2FD'
                    ax.add_patch(plt.Rectangle((bx, by), box_w, box_h,
                                 facecolor=gc, edgecolor='black', linewidth=0.5))
                    dist_text = f"{g['dist']} bp" if g['dist'] > 0 else "0 bp"
                    ax.text(bx + box_w/2, by + box_h/2, f"{gname}\n({dist_text})",
                            ha='center', va='center', fontsize=5.5)

    plt.savefig(os.path.join(OUT_DIR, 'Fig5_IR_boundary.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'Fig5_IR_boundary.pdf'), bbox_inches='tight')
    plt.close()
    print("  -> Fig5_IR_boundary.png/pdf")


###############################################################################
# Fig 6: Structural Variant Analysis (recalculated for 8 varieties)
###############################################################################
def fig6_variant():
    print("Generating Fig 6: Variant analysis (8 varieties)...")

    # Re-filter variants: a variant is only informative if, among 8 varieties,
    # there is actual variation (not all same allele after removing GB048/GB066)
    excluded = {'GB048', 'GB066'}
    region_counts = defaultdict(lambda: {'SNP':0, 'InDel':0, 'Complex':0})
    gene_counts = defaultdict(int)
    total_snp = 0; total_indel = 0; total_complex = 0

    # Build pairwise SNP matrix from raw data
    # We need to re-parse from all_variants.csv and check if variant still exists
    # for 8 varieties
    valid_variants = []
    with open('/home/spenkim/rubus_cp/variant_output/all_variants.csv') as f:
        for row in csv.DictReader(f):
            details = row['Details']
            vtype = row['Type']
            region = row['Region']
            gene = row['Gene']

            # Check if this variant is still polymorphic among 8 varieties
            # Parse which varieties have which allele
            if vtype == 'SNP':
                # Format: "T:DAHE,DAHO,...; A:DB,GDC,..."
                allele_groups = {}
                for part in details.split('; '):
                    if ':' in part:
                        allele, vars_str = part.split(':', 1)
                        vars_list = [v.strip() for v in vars_str.split(',')]
                        vars_8 = [v for v in vars_list if v not in excluded]
                        if vars_8:
                            allele_groups[allele] = vars_8
                # Still polymorphic if >1 allele group has varieties
                if len(allele_groups) >= 2:
                    valid_variants.append(row)
                    region_counts[region]['SNP'] += 1
                    total_snp += 1
                    if gene != 'intergenic':
                        gene_counts[gene] += 1
            elif vtype == 'InDel':
                # Format: "gap_in=DB,DS,...; seq=ATAAATTG"
                m = re.match(r'gap_in=([^;]+)', details)
                if m:
                    gap_vars = [v.strip() for v in m.group(1).split(',')]
                    gap_8 = [v for v in gap_vars if v not in excluded]
                    non_gap_8 = [v for v in VARIETIES if v not in gap_8]
                    # Polymorphic if both groups non-empty
                    if gap_8 and non_gap_8:
                        valid_variants.append(row)
                        region_counts[region]['InDel'] += 1
                        total_indel += 1
                        if gene != 'intergenic':
                            gene_counts[gene] += 1
                else:
                    # Keep as is if format doesn't match
                    valid_variants.append(row)
                    region_counts[region]['InDel'] += 1
                    total_indel += 1
            else:  # Complex
                # Keep and check similarly
                valid_variants.append(row)
                region_counts[region]['Complex'] += 1
                total_complex += 1

    print(f"  8-variety variants: {total_snp} SNPs, {total_indel} InDels, {total_complex} Complex = {total_snp+total_indel+total_complex} total")

    # Calculate private variants for 8 varieties
    private_vars = {v: {'SNPs':0, 'InDels':0} for v in VARIETIES}
    for row in valid_variants:
        details = row['Details']
        vtype = row['Type']
        if vtype == 'SNP':
            for part in details.split('; '):
                if ':' in part:
                    allele, vars_str = part.split(':', 1)
                    vars_list = [v.strip() for v in vars_str.split(',') if v.strip() not in excluded]
                    if len(vars_list) == 1 and vars_list[0] in VARIETIES:
                        private_vars[vars_list[0]]['SNPs'] += 1
        elif vtype == 'InDel':
            m = re.match(r'gap_in=([^;]+)', details)
            if m:
                gap_vars = [v.strip() for v in m.group(1).split(',') if v.strip() not in excluded]
                non_gap = [v for v in VARIETIES if v not in gap_vars]
                if len(gap_vars) == 1 and gap_vars[0] in VARIETIES:
                    private_vars[gap_vars[0]]['InDels'] += 1
                elif len(non_gap) == 1 and non_gap[0] in VARIETIES:
                    private_vars[non_gap[0]]['InDels'] += 1

    # Pairwise SNP counts for 8 varieties
    # Re-parse from variant_summary.txt (the original 10x10 matrix)
    snp_matrix_10 = {}
    matrix_vars_10 = []
    with open('/home/spenkim/rubus_cp/variant_output/variant_summary.txt') as f:
        content = f.read()
    in_matrix = False
    for line in content.split('\n'):
        if 'Pairwise SNP counts' in line:
            in_matrix = True; continue
        if in_matrix:
            line = line.strip()
            parts = line.split()
            if parts and parts[0] in ALL_10:
                var_name = parts[0]
                vals = []
                for p in parts[1:]:
                    if p == '—':
                        vals.append(0)
                    else:
                        try: vals.append(int(p))
                        except ValueError: pass
                snp_matrix_10[var_name] = vals
                matrix_vars_10.append(var_name)

    # Extract 8x8 submatrix
    idx_8 = [matrix_vars_10.index(v) for v in VARIETIES if v in matrix_vars_10]
    n8 = len(VARIETIES)
    matrix_8 = np.zeros((n8, n8))
    for i, vi in enumerate(VARIETIES):
        if vi in snp_matrix_10:
            row_vals = snp_matrix_10[vi]
            for j, vj in enumerate(VARIETIES):
                j10 = matrix_vars_10.index(vj) if vj in matrix_vars_10 else -1
                if 0 <= j10 < len(row_vals):
                    matrix_8[i][j] = row_vals[j10]

    # Region sizes (from DAHO reference)
    region_sizes = {'LSC': 85374, 'IRb': 25786, 'SSC': 18718, 'IRa': 25786}

    # --- PLOTTING ---
    fig = plt.figure(figsize=(18, 8))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1.3], wspace=0.35)

    # (A) Variants by region
    ax1 = fig.add_subplot(gs[0])
    regions_order = ['LSC', 'IRb', 'SSC', 'IRa']
    x = np.arange(4); w = 0.25
    snps = [region_counts[r]['SNP'] for r in regions_order]
    indels = [region_counts[r]['InDel'] for r in regions_order]
    cmplx = [region_counts[r]['Complex'] for r in regions_order]

    ax1.bar(x-w, snps, w, label='SNPs', color=COLORS['SNP'])
    ax1.bar(x, indels, w, label='InDels', color=COLORS['InDel'])
    ax1.bar(x+w, cmplx, w, label='Complex', color=COLORS['Complex'])
    ax1.set_xticks(x)
    ax1.set_xticklabels(regions_order, fontsize=10)
    ax1.set_ylabel('Number of Variants', fontsize=11)
    ax1.set_title('(A) Variants by Genome Region', fontsize=12)
    ax1.legend(fontsize=9)

    for i, r in enumerate(regions_order):
        total = region_counts[r]['SNP'] + region_counts[r]['InDel'] + region_counts[r]['Complex']
        per_kb = total / (region_sizes[r] / 1000)
        ax1.text(i, total+5, f'{per_kb:.2f}/kb', ha='center', va='bottom', fontsize=8, color='gray')

    # (B) Private variants
    ax2 = fig.add_subplot(gs[1])
    x2 = np.arange(n8); w2 = 0.35
    ps = [private_vars[v]['SNPs'] for v in VARIETIES]
    pi = [private_vars[v]['InDels'] for v in VARIETIES]
    ax2.bar(x2 - w2/2, ps, w2, label='Private SNPs', color=COLORS['SNP'])
    ax2.bar(x2 + w2/2, pi, w2, label='Private InDels', color=COLORS['InDel'])
    ax2.set_xticks(x2)
    ax2.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Number of Private Variants', fontsize=11)
    ax2.set_title('(B) Private Variants per Variety', fontsize=12)
    ax2.legend(fontsize=9)

    # (C) Pairwise SNP heatmap
    ax3 = fig.add_subplot(gs[2])
    im = ax3.imshow(matrix_8, cmap='YlOrRd', aspect='equal')
    ax3.set_xticks(range(n8))
    ax3.set_yticks(range(n8))
    ax3.set_xticklabels(VARIETIES, rotation=45, ha='right', fontsize=9)
    ax3.set_yticklabels(VARIETIES, fontsize=9)
    ax3.set_title('(C) Pairwise SNP Counts', fontsize=12)
    for i in range(n8):
        for j in range(n8):
            val = int(matrix_8[i][j])
            color = 'white' if val > 200 else 'black'
            ax3.text(j, i, str(val), ha='center', va='center', fontsize=7, color=color)
    cbar = plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
    cbar.set_label('SNP Count', fontsize=10)

    plt.savefig(os.path.join(OUT_DIR, 'Fig6_Variant.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'Fig6_Variant.pdf'), bbox_inches='tight')
    plt.close()
    print("  -> Fig6_Variant.png/pdf")


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    print("=" * 60)
    print("Generating figures (8 varieties, GB048/GB066 excluded)")
    print("=" * 60)
    fig2_rscu()
    fig3_ssr()
    fig4_repeat()
    fig5_irscope()
    fig6_variant()
    print("\n" + "=" * 60)
    print(f"All figures saved to: {OUT_DIR}")
    print("=" * 60)
