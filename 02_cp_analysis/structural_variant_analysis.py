#!/usr/bin/env python3
"""
structural_variant_analysis.py — SNP & InDel detection among 10 Rubus cp genomes

Uses MAFFT alignment of the 10 varieties to detect:
1. SNPs (single nucleotide polymorphisms)
2. InDels (insertions/deletions)
3. Variant distribution across genome regions (LSC, IRb, SSC, IRa)
4. Per-gene variant annotation
"""

import os
import csv
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

ALIGNMENT = f"{BASE}/standardized/all10_std_mafft.fasta"
REFERENCE_GB = f"{BASE}/bankit_gb/DAHO.gb"  # Use DAHO as reference
OUT_DIR = f"{BASE}/variant_output"

VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]


def load_alignment():
    """Load MAFFT alignment."""
    records = {}
    for rec in SeqIO.parse(ALIGNMENT, "fasta"):
        records[rec.id] = str(rec.seq).upper()
    return records


def detect_variants(records):
    """Detect SNPs and InDels from alignment.

    Returns list of variant dicts.
    """
    seqnames = list(records.keys())
    aln_len = len(records[seqnames[0]])
    variants = []

    i = 0
    while i < aln_len:
        bases = {name: records[name][i] for name in seqnames}
        unique = set(bases.values())

        # All identical
        if len(unique) == 1:
            i += 1
            continue

        # Check if it's a gap-involved position (InDel) or pure SNP
        has_gap = "-" in unique
        non_gap = unique - {"-"}

        if has_gap and len(non_gap) <= 1:
            # Pure InDel — find extent of contiguous gap block
            indel_start = i
            # Collect which samples have gaps
            gap_samples = {name for name in seqnames if bases[name] == "-"}
            non_gap_samples = set(seqnames) - gap_samples

            # Extend while the same gap pattern continues
            j = i + 1
            while j < aln_len:
                next_bases = {name: records[name][j] for name in seqnames}
                next_gap = {name for name in seqnames if next_bases[name] == "-"}
                if next_gap == gap_samples:
                    j += 1
                else:
                    break

            indel_len = j - i
            # Get the inserted sequence from non-gap samples
            ref_sample = non_gap_samples.pop()
            inserted_seq = records[ref_sample][i:j].replace("-", "")

            variants.append({
                "type": "InDel",
                "aln_start": i,
                "aln_end": j - 1,
                "length": indel_len,
                "gap_samples": sorted(gap_samples),
                "non_gap_samples": sorted(set(seqnames) - gap_samples),
                "inserted_seq": inserted_seq[:50],  # truncate long seqs
            })
            i = j
            continue

        elif has_gap and len(non_gap) > 1:
            # Mixed: gap + multiple alleles — treat as complex variant
            # Find extent
            indel_start = i
            j = i + 1
            while j < aln_len:
                next_bases = {name: records[name][j] for name in seqnames}
                next_unique = set(next_bases.values())
                if "-" in next_unique:
                    j += 1
                else:
                    break

            alleles = defaultdict(list)
            for name in seqnames:
                allele = records[name][i:j].replace("-", "")
                alleles[allele].append(name)

            variants.append({
                "type": "Complex",
                "aln_start": i,
                "aln_end": j - 1,
                "length": j - i,
                "alleles": dict(alleles),
            })
            i = j
            continue

        else:
            # Pure SNP (no gaps)
            alleles = defaultdict(list)
            for name in seqnames:
                alleles[bases[name]].append(name)

            variants.append({
                "type": "SNP",
                "aln_start": i,
                "aln_end": i,
                "length": 1,
                "alleles": dict(alleles),
            })
            i += 1
            continue

    return variants


def aln_to_ref_pos(records, ref_name="DAHO"):
    """Build alignment position → reference genome position mapping."""
    ref_seq = records[ref_name]
    mapping = {}
    ref_pos = 0
    for aln_pos in range(len(ref_seq)):
        if ref_seq[aln_pos] != "-":
            mapping[aln_pos] = ref_pos
            ref_pos += 1
        else:
            mapping[aln_pos] = None  # gap in reference
    return mapping


def get_genome_regions():
    """Get LSC/IRb/SSC/IRa boundaries from reference GenBank."""
    record = SeqIO.read(REFERENCE_GB, "genbank")
    regions = {}
    for feat in record.features:
        if feat.type == "misc_feature":
            note = feat.qualifiers.get("note", [""])[0].lower()
            if "large single copy" in note:
                regions["LSC"] = (int(feat.location.start), int(feat.location.end))
            elif "small single copy" in note:
                regions["SSC"] = (int(feat.location.start), int(feat.location.end))
        elif feat.type == "repeat_region":
            note = feat.qualifiers.get("note", [""])[0].lower()
            if "irb" in note or "inverted repeat b" in note:
                regions["IRb"] = (int(feat.location.start), int(feat.location.end))
            elif "ira" in note or "inverted repeat a" in note:
                regions["IRa"] = (int(feat.location.start), int(feat.location.end))

    # Fallback
    if len(regions) < 4:
        rr = []
        for feat in record.features:
            if feat.type == "repeat_region":
                rr.append((int(feat.location.start), int(feat.location.end)))
        if len(rr) == 2:
            rr.sort()
            regions["IRb"] = rr[0]
            regions["IRa"] = rr[1]
        mf = []
        for feat in record.features:
            if feat.type == "misc_feature":
                mf.append((int(feat.location.start), int(feat.location.end)))
        if len(mf) == 2:
            mf.sort()
            regions["LSC"] = mf[0]
            regions["SSC"] = mf[1]

    return regions


def get_gene_map():
    """Build position → gene name mapping from reference GenBank.

    For trans-spliced genes (like rps12), use individual exon/CDS parts
    instead of the full gene span to avoid inflated variant counts.
    """
    record = SeqIO.read(REFERENCE_GB, "genbank")
    gene_intervals = []
    seen = set()

    # First pass: collect CDS/tRNA/rRNA features (actual coding regions)
    coding_intervals = set()
    for feat in record.features:
        if feat.type not in ("CDS", "tRNA", "rRNA"):
            continue
        name = feat.qualifiers.get("gene", ["?"])[0]
        # Use individual parts for compound locations (trans-spliced genes)
        for part in feat.location.parts:
            key = f"{name}_{int(part.start)}_{int(part.end)}"
            if key not in seen:
                seen.add(key)
                gene_intervals.append((int(part.start), int(part.end), name))
                coding_intervals.add(key)

    # Second pass: add gene features that don't have CDS/tRNA/rRNA
    # (e.g., pseudogenes with only gene feature)
    gene_names_with_coding = {name for _, _, name in gene_intervals}
    for feat in record.features:
        if feat.type != "gene":
            continue
        name = feat.qualifiers.get("gene", ["?"])[0]
        # For genes that already have coding features, skip (they're covered by parts above)
        # For simple genes without compound locations that aren't trans-spliced, add them
        if name in gene_names_with_coding:
            continue
        start = int(feat.location.start)
        end = int(feat.location.end)
        key = f"{name}_{start}_{end}"
        if key not in seen:
            seen.add(key)
            gene_intervals.append((start, end, name))

    return sorted(gene_intervals)


def pos_to_region(pos, regions):
    """Map a genome position to a region name."""
    for rname, (rstart, rend) in regions.items():
        if rstart <= pos < rend:
            return rname
    return "unknown"


def pos_to_gene(pos, gene_intervals):
    """Map a genome position to a gene name (or 'intergenic')."""
    for gstart, gend, gname in gene_intervals:
        if gstart <= pos < gend:
            return gname
    return "intergenic"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Structural Variant Analysis")
    print("=" * 40)

    # Load data
    print("Loading alignment...")
    records = load_alignment()
    print(f"  {len(records)} sequences, {len(next(iter(records.values()))):,} aligned positions")

    # Detect variants
    print("Detecting variants...")
    variants = detect_variants(records)
    print(f"  Total variants: {len(variants)}")

    snps = [v for v in variants if v["type"] == "SNP"]
    indels = [v for v in variants if v["type"] == "InDel"]
    complex_vars = [v for v in variants if v["type"] == "Complex"]
    print(f"  SNPs: {len(snps)}")
    print(f"  InDels: {len(indels)}")
    print(f"  Complex: {len(complex_vars)}")

    # Map to reference positions
    print("Mapping to reference genome positions...")
    pos_map = aln_to_ref_pos(records, "DAHO")
    regions = get_genome_regions()
    gene_intervals = get_gene_map()

    # Annotate variants with region and gene
    for v in variants:
        ref_pos = pos_map.get(v["aln_start"])
        if ref_pos is not None:
            v["ref_pos"] = ref_pos
            v["region"] = pos_to_region(ref_pos, regions)
            v["gene"] = pos_to_gene(ref_pos, gene_intervals)
        else:
            # Position is a gap in reference (insertion relative to DAHO)
            # Find nearest non-gap position
            for offset in range(1, 100):
                if pos_map.get(v["aln_start"] - offset) is not None:
                    ref_pos = pos_map[v["aln_start"] - offset]
                    break
            v["ref_pos"] = ref_pos if ref_pos is not None else -1
            v["region"] = pos_to_region(v["ref_pos"], regions) if v["ref_pos"] >= 0 else "unknown"
            v["gene"] = pos_to_gene(v["ref_pos"], gene_intervals) if v["ref_pos"] >= 0 else "unknown"

    # --- Output files ---

    # 1. Full variant list CSV
    csv_path = os.path.join(OUT_DIR, "all_variants.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "Ref_pos", "Aln_start", "Aln_end", "Length",
                         "Region", "Gene", "Details"])
        for v in variants:
            if v["type"] == "SNP":
                details = "; ".join(f"{allele}:{','.join(samples)}"
                                   for allele, samples in v["alleles"].items())
            elif v["type"] == "InDel":
                details = (f"gap_in={','.join(v['gap_samples'])}; "
                          f"seq={v['inserted_seq']}")
            else:
                details = "; ".join(f"{allele or 'DEL'}:{','.join(samples)}"
                                   for allele, samples in v["alleles"].items())
            writer.writerow([v["type"], v.get("ref_pos", ""), v["aln_start"],
                           v["aln_end"], v["length"], v["region"], v["gene"], details])
    print(f"\nAll variants: {csv_path}")

    # 2. Summary by region
    summary_path = os.path.join(OUT_DIR, "variant_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Structural Variant Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Alignment: {ALIGNMENT}\n")
        f.write(f"Reference: DAHO\n")
        f.write(f"Varieties: {', '.join(VARIETIES)}\n\n")

        f.write(f"Total variants: {len(variants)}\n")
        f.write(f"  SNPs: {len(snps)}\n")
        f.write(f"  InDels: {len(indels)}\n")
        f.write(f"  Complex: {len(complex_vars)}\n\n")

        # By region
        f.write("Variants by genome region:\n")
        f.write(f"{'Region':<10}{'SNPs':>8}{'InDels':>8}{'Complex':>8}{'Total':>8}\n")
        f.write("-" * 42 + "\n")
        for rname in ["LSC", "IRb", "SSC", "IRa", "unknown"]:
            n_snp = sum(1 for v in snps if v.get("region") == rname)
            n_indel = sum(1 for v in indels if v.get("region") == rname)
            n_complex = sum(1 for v in complex_vars if v.get("region") == rname)
            total = n_snp + n_indel + n_complex
            if total > 0:
                f.write(f"{rname:<10}{n_snp:>8}{n_indel:>8}{n_complex:>8}{total:>8}\n")
        f.write("\n")

        # By region normalized (per kb)
        f.write("Variant density (per kb):\n")
        region_sizes = {r: (e - s) for r, (s, e) in regions.items()}
        f.write(f"{'Region':<10}{'Size_bp':>10}{'Variants':>10}{'Per_kb':>10}\n")
        f.write("-" * 40 + "\n")
        for rname in ["LSC", "IRb", "SSC", "IRa"]:
            n = sum(1 for v in variants if v.get("region") == rname)
            size = region_sizes.get(rname, 1)
            per_kb = n / (size / 1000)
            f.write(f"{rname:<10}{size:>10,}{n:>10}{per_kb:>10.2f}\n")
        f.write("\n")

        # Genic vs intergenic
        n_genic = sum(1 for v in variants if v.get("gene", "intergenic") != "intergenic")
        n_intergenic = sum(1 for v in variants if v.get("gene", "intergenic") == "intergenic")
        f.write(f"Genic variants: {n_genic}\n")
        f.write(f"Intergenic variants: {n_intergenic}\n\n")

        # Top variant genes
        gene_counts = Counter(v.get("gene", "intergenic") for v in variants if v.get("gene") != "intergenic")
        f.write("Top 20 genes with most variants:\n")
        f.write(f"{'Gene':<20}{'Variants':>10}\n")
        f.write("-" * 30 + "\n")
        for gene, count in gene_counts.most_common(20):
            f.write(f"{gene:<20}{count:>10}\n")
        f.write("\n")

        # InDel size distribution
        f.write("InDel size distribution:\n")
        indel_sizes = [v["length"] for v in indels]
        if indel_sizes:
            size_bins = [(1, 1), (2, 5), (6, 10), (11, 50), (51, 100), (101, 500)]
            f.write(f"{'Size_range':<15}{'Count':>8}\n")
            f.write("-" * 23 + "\n")
            for lo, hi in size_bins:
                n = sum(1 for s in indel_sizes if lo <= s <= hi)
                if n > 0:
                    f.write(f"{lo}-{hi} bp{'':<6}{n:>8}\n")
            over500 = sum(1 for s in indel_sizes if s > 500)
            if over500:
                f.write(f">500 bp{'':<8}{over500:>8}\n")
            f.write(f"\nMin InDel: {min(indel_sizes)} bp\n")
            f.write(f"Max InDel: {max(indel_sizes)} bp\n")
            f.write(f"Mean InDel: {sum(indel_sizes)/len(indel_sizes):.1f} bp\n")

        # Pairwise SNP counts
        f.write("\n\nPairwise SNP counts:\n")
        names = sorted(records.keys())
        f.write(f"{'':>8}" + "".join(f"{n:>8}" for n in names) + "\n")
        for n1 in names:
            row = f"{n1:>8}"
            for n2 in names:
                if n1 == n2:
                    row += f"{'—':>8}"
                else:
                    count = 0
                    for v in snps:
                        alleles = v.get("alleles", {})
                        a1 = a2 = None
                        for allele, samples in alleles.items():
                            if n1 in samples:
                                a1 = allele
                            if n2 in samples:
                                a2 = allele
                        if a1 and a2 and a1 != a2:
                            count += 1
                    row += f"{count:>8}"
            f.write(row + "\n")

    print(f"Summary: {summary_path}")

    # 3. Per-variety variant count
    var_count_path = os.path.join(OUT_DIR, "per_variety_variants.csv")
    with open(var_count_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Variety", "Private_SNPs", "Private_InDels"])

        for name in VARIETIES:
            # Private SNP: this variety has a unique allele
            priv_snp = 0
            for v in snps:
                for allele, samples in v["alleles"].items():
                    if samples == [name]:
                        priv_snp += 1

            # Private InDel: gap pattern unique to this variety
            priv_indel = 0
            for v in indels:
                if v["gap_samples"] == [name] or v["non_gap_samples"] == [name]:
                    priv_indel += 1

            writer.writerow([name, priv_snp, priv_indel])

    print(f"Per-variety: {var_count_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
