#!/usr/bin/env python3
"""
rscu_analysis.py — Codon Usage & RSCU analysis for Rubus cp genomes

Extracts all protein-coding CDS from GenBank files,
calculates codon frequencies and RSCU values.
Outputs per-variety and summary tables (TSV + CSV).
"""

import os
import csv
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import unambiguous_dna_by_id

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]
GB_DIR = f"{BASE}/bankit_gb"
OUT_DIR = f"{BASE}/rscu_output"

# Plastid genetic code (table 11 = bacterial/plant plastid)
CODON_TABLE = unambiguous_dna_by_id[11]

# Standard 64 codons
BASES = "TCAG"
ALL_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]

# Build amino acid -> synonymous codons mapping
AA_TO_CODONS = defaultdict(list)
for codon in ALL_CODONS:
    if codon in CODON_TABLE.stop_codons:
        AA_TO_CODONS["*"].append(codon)
    else:
        aa = CODON_TABLE.forward_table[codon]
        AA_TO_CODONS[aa].append(codon)

# Codon to amino acid
CODON_TO_AA = {}
for codon in ALL_CODONS:
    if codon in CODON_TABLE.stop_codons:
        CODON_TO_AA[codon] = "*"
    else:
        CODON_TO_AA[codon] = CODON_TABLE.forward_table[codon]


def extract_cds_sequences(gb_path):
    """Extract CDS nucleotide sequences from a GenBank file.

    Returns list of (gene_name, dna_sequence) tuples.
    Skips pseudogenes and IR duplicate CDS (keeps only first occurrence).
    """
    record = SeqIO.read(gb_path, "genbank")
    seen_genes = set()
    cds_list = []

    for feat in record.features:
        if feat.type != "CDS":
            continue

        # Skip pseudogenes
        if "pseudo" in feat.qualifiers:
            continue

        gene = feat.qualifiers.get("gene", ["unknown"])[0]

        # Skip IR duplicates (keep first copy only to avoid double-counting)
        if gene in seen_genes:
            continue
        seen_genes.add(gene)

        # Extract sequence
        seq = feat.location.extract(record.seq)
        dna = str(seq).upper()

        # Validate: length must be multiple of 3
        if len(dna) % 3 != 0:
            print(f"  WARNING: {gene} length {len(dna)} not divisible by 3, skipping")
            continue

        cds_list.append((gene, dna))

    return cds_list


def count_codons(cds_list):
    """Count codon usage across all CDS sequences.

    Excludes stop codons and start codons (first codon of each CDS).
    """
    counts = Counter()
    total_codons = 0

    for gene, dna in cds_list:
        # Process codons (skip first = start codon, skip last = stop codon)
        for i in range(3, len(dna) - 3, 3):
            codon = dna[i:i+3]
            if 'N' in codon:
                continue
            counts[codon] += 1
            total_codons += 1

    return counts, total_codons


def calculate_rscu(codon_counts):
    """Calculate RSCU (Relative Synonymous Codon Usage).

    RSCU = observed frequency / expected frequency
    Expected = total_for_aa / number_of_synonymous_codons
    RSCU = (count * num_synonymous) / total_for_aa

    If total_for_aa = 0, RSCU = 0.
    """
    rscu = {}
    for aa, codons in AA_TO_CODONS.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        n_syn = len(codons)
        for c in codons:
            if total == 0:
                rscu[c] = 0.0
            else:
                rscu[c] = (codon_counts.get(c, 0) * n_syn) / total
    return rscu


def analyze_variety(variety):
    """Run full codon analysis for one variety."""
    gb_path = os.path.join(GB_DIR, f"{variety}.gb")
    cds_list = extract_cds_sequences(gb_path)
    codon_counts, total_codons = count_codons(cds_list)
    rscu = calculate_rscu(codon_counts)
    return cds_list, codon_counts, total_codons, rscu


def write_rscu_table(all_results):
    """Write comprehensive RSCU table (all varieties)."""
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Per-variety RSCU CSV
    rscu_path = os.path.join(OUT_DIR, "rscu_all_varieties.csv")
    with open(rscu_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Codon", "AminoAcid"] + VARIETIES + ["Mean_RSCU"]
        writer.writerow(header)

        for codon in ALL_CODONS:
            aa = CODON_TO_AA[codon]
            row = [codon, aa]
            values = []
            for var in VARIETIES:
                v = all_results[var]["rscu"].get(codon, 0.0)
                row.append(f"{v:.3f}")
                values.append(v)
            row.append(f"{sum(values)/len(values):.3f}")
            writer.writerow(row)

    print(f"  RSCU table: {rscu_path}")

    # 2. Per-variety codon count CSV
    count_path = os.path.join(OUT_DIR, "codon_counts_all_varieties.csv")
    with open(count_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Codon", "AminoAcid"] + VARIETIES
        writer.writerow(header)

        for codon in ALL_CODONS:
            aa = CODON_TO_AA[codon]
            row = [codon, aa]
            for var in VARIETIES:
                row.append(all_results[var]["counts"].get(codon, 0))
            writer.writerow(row)

    print(f"  Codon counts: {count_path}")

    # 3. Summary statistics
    summary_path = os.path.join(OUT_DIR, "rscu_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Codon Usage Analysis Summary\n")
        f.write("=" * 60 + "\n\n")

        for var in VARIETIES:
            r = all_results[var]
            f.write(f"{var}:\n")
            f.write(f"  CDS analyzed: {r['n_cds']}\n")
            f.write(f"  Total codons (excl. start/stop): {r['total_codons']}\n")

            # Count preferred codons (RSCU > 1)
            preferred = sum(1 for c in ALL_CODONS
                          if CODON_TO_AA[c] != "*" and r["rscu"].get(c, 0) > 1.0)
            f.write(f"  Preferred codons (RSCU > 1): {preferred}\n")

            # A/T ending bias
            at_end = sum(r["counts"].get(c, 0) for c in ALL_CODONS
                        if c[2] in "AT" and CODON_TO_AA[c] != "*")
            gc_end = sum(r["counts"].get(c, 0) for c in ALL_CODONS
                        if c[2] in "GC" and CODON_TO_AA[c] != "*")
            total = at_end + gc_end
            if total > 0:
                f.write(f"  A/T-ending codons: {at_end} ({100*at_end/total:.1f}%)\n")
                f.write(f"  G/C-ending codons: {gc_end} ({100*gc_end/total:.1f}%)\n")
            f.write("\n")

        # Mean RSCU table (publication format)
        f.write("\n" + "=" * 60 + "\n")
        f.write("Mean RSCU across all 10 varieties (publication format)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"{'Codon':<8}{'AA':<6}{'Count':>8}{'RSCU':>8}\n")
        f.write("-" * 30 + "\n")

        # Use DAHO as representative for counts (all very similar)
        ref = all_results["DAHO"]
        for aa in sorted(set(CODON_TO_AA.values())):
            if aa == "*":
                continue
            codons = AA_TO_CODONS[aa]
            for c in sorted(codons):
                mean_rscu = sum(all_results[v]["rscu"].get(c, 0) for v in VARIETIES) / len(VARIETIES)
                mean_count = sum(all_results[v]["counts"].get(c, 0) for v in VARIETIES) / len(VARIETIES)
                f.write(f"{c:<8}{aa:<6}{mean_count:>8.0f}{mean_rscu:>8.3f}\n")
            f.write("\n")

    print(f"  Summary: {summary_path}")

    # 4. DAHO representative detailed table (for paper Table)
    paper_path = os.path.join(OUT_DIR, "rscu_paper_table.csv")
    with open(paper_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AminoAcid", "Codon", "Number", "RSCU"])

        ref_counts = all_results["DAHO"]["counts"]
        ref_rscu = all_results["DAHO"]["rscu"]

        for aa in sorted(set(CODON_TO_AA.values())):
            if aa == "*":
                continue
            codons = AA_TO_CODONS[aa]
            for c in sorted(codons):
                writer.writerow([
                    aa, c,
                    ref_counts.get(c, 0),
                    f"{ref_rscu.get(c, 0.0):.3f}"
                ])

    print(f"  Paper table (DAHO): {paper_path}")


def main():
    print("Codon Usage & RSCU Analysis")
    print("=" * 40)

    all_results = {}
    for var in VARIETIES:
        print(f"\n{var}:")
        cds_list, counts, total, rscu = analyze_variety(var)
        all_results[var] = {
            "cds_list": cds_list,
            "n_cds": len(cds_list),
            "counts": counts,
            "total_codons": total,
            "rscu": rscu,
        }
        print(f"  {len(cds_list)} CDS, {total} codons analyzed")

    print("\nWriting output files...")
    write_rscu_table(all_results)
    print("\nDone.")


if __name__ == "__main__":
    main()
