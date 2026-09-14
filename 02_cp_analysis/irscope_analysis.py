#!/usr/bin/env python3
"""
irscope_analysis.py — IR boundary gene analysis for Rubus cp genomes

Identifies genes at the four IR junctions (JLB, JSB, JSA, JLA) for each variety.
Produces a tabular summary for IRscope-style figures.

Junctions:
  JLB = LSC/IRb boundary
  JSB = IRb/SSC boundary
  JSA = SSC/IRa boundary
  JLA = IRa/LSC boundary (= genome end/start for circular)
"""

import os
import csv
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]
GB_DIR = f"{BASE}/bankit_gb"
OUT_DIR = f"{BASE}/irscope_output"


def get_regions(record):
    """Extract LSC, IRb, SSC, IRa boundaries from misc_feature/repeat_region."""
    regions = {}
    for feat in record.features:
        if feat.type == "misc_feature":
            note = feat.qualifiers.get("note", [""])[0].lower()
            if "large single copy" in note or "lsc" in note:
                regions["LSC"] = (int(feat.location.start), int(feat.location.end))
            elif "small single copy" in note or "ssc" in note:
                regions["SSC"] = (int(feat.location.start), int(feat.location.end))
        elif feat.type == "repeat_region":
            note = feat.qualifiers.get("note", [""])[0].lower()
            rpt = feat.qualifiers.get("rpt_type", [""])[0].lower()
            if "irb" in note or "inverted repeat b" in note or (rpt == "inverted" and "b" in note):
                regions["IRb"] = (int(feat.location.start), int(feat.location.end))
            elif "ira" in note or "inverted repeat a" in note or (rpt == "inverted" and "a" in note):
                regions["IRa"] = (int(feat.location.start), int(feat.location.end))

    # If labels didn't work, try by position order
    if len(regions) < 4:
        # Fallback: find repeat_regions and assign by position
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


def get_gene_features(record):
    """Get all gene features with their positions."""
    genes = []
    for feat in record.features:
        if feat.type != "gene":
            continue
        gene_name = feat.qualifiers.get("gene", ["?"])[0]
        is_pseudo = "pseudo" in feat.qualifiers

        # Get overall span
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand

        genes.append({
            "name": gene_name,
            "start": start,
            "end": end,
            "strand": strand,
            "pseudo": is_pseudo,
            "length": end - start,
        })
    return genes


def find_junction_genes(genes, regions, genome_len):
    """Find genes at each of the four IR junctions."""
    lsc_start, lsc_end = regions["LSC"]
    irb_start, irb_end = regions["IRb"]
    ssc_start, ssc_end = regions["SSC"]
    ira_start, ira_end = regions["IRa"]

    # Junction positions (0-based)
    jlb = irb_start   # LSC/IRb
    jsb = irb_end      # IRb/SSC
    jsa = ira_start    # SSC/IRa
    jla = ira_end      # IRa/LSC (wraps to position 0)

    results = {}

    # --- JLB (LSC/IRb) ---
    # Find last gene ending in LSC near junction & first gene starting in IRb
    jlb_genes = []
    for g in genes:
        # Gene spanning the junction
        if g["start"] < jlb < g["end"]:
            dist_in_lsc = jlb - g["start"]
            dist_in_irb = g["end"] - jlb
            jlb_genes.append({**g, "position": "spanning",
                             "bp_in_LSC": dist_in_lsc, "bp_in_IRb": dist_in_irb,
                             "dist_to_junction": 0})
        # Gene entirely in LSC, close to junction
        elif g["end"] <= jlb and (jlb - g["end"]) < 5000:
            jlb_genes.append({**g, "position": "LSC",
                             "dist_to_junction": jlb - g["end"]})
        # Gene entirely in IRb, close to junction
        elif g["start"] >= jlb and (g["start"] - jlb) < 5000:
            jlb_genes.append({**g, "position": "IRb",
                             "dist_to_junction": g["start"] - jlb})

    results["JLB"] = sorted(jlb_genes, key=lambda x: abs(x.get("dist_to_junction", 0)))[:4]

    # --- JSB (IRb/SSC) ---
    jsb_genes = []
    for g in genes:
        if g["start"] < jsb < g["end"]:
            dist_in_irb = jsb - g["start"]
            dist_in_ssc = g["end"] - jsb
            jsb_genes.append({**g, "position": "spanning",
                             "bp_in_IRb": dist_in_irb, "bp_in_SSC": dist_in_ssc,
                             "dist_to_junction": 0})
        elif g["end"] <= jsb and (jsb - g["end"]) < 5000:
            jsb_genes.append({**g, "position": "IRb",
                             "dist_to_junction": jsb - g["end"]})
        elif g["start"] >= jsb and (g["start"] - jsb) < 5000:
            jsb_genes.append({**g, "position": "SSC",
                             "dist_to_junction": g["start"] - jsb})

    results["JSB"] = sorted(jsb_genes, key=lambda x: abs(x.get("dist_to_junction", 0)))[:4]

    # --- JSA (SSC/IRa) ---
    jsa_genes = []
    for g in genes:
        if g["start"] < jsa < g["end"]:
            dist_in_ssc = jsa - g["start"]
            dist_in_ira = g["end"] - jsa
            jsa_genes.append({**g, "position": "spanning",
                             "bp_in_SSC": dist_in_ssc, "bp_in_IRa": dist_in_ira,
                             "dist_to_junction": 0})
        elif g["end"] <= jsa and (jsa - g["end"]) < 5000:
            jsa_genes.append({**g, "position": "SSC",
                             "dist_to_junction": jsa - g["end"]})
        elif g["start"] >= jsa and (g["start"] - jsa) < 5000:
            jsa_genes.append({**g, "position": "IRa",
                             "dist_to_junction": g["start"] - jsa})

    results["JSA"] = sorted(jsa_genes, key=lambda x: abs(x.get("dist_to_junction", 0)))[:4]

    # --- JLA (IRa/LSC) — wraps around genome start ---
    jla_genes = []
    for g in genes:
        # Gene at end of IRa (close to genome end)
        if g["end"] <= genome_len and (genome_len - g["end"]) < 5000 and g["start"] >= ira_start:
            jla_genes.append({**g, "position": "IRa",
                             "dist_to_junction": genome_len - g["end"]})
        # Gene at start of genome (LSC start)
        elif g["start"] < 5000 and g["end"] < lsc_end:
            jla_genes.append({**g, "position": "LSC",
                             "dist_to_junction": g["start"]})
        # Gene spanning genome end/start
        if g["start"] >= ira_start and g["end"] > genome_len - 100:
            # Check if it extends past the genome boundary
            pass

    results["JLA"] = sorted(jla_genes, key=lambda x: abs(x.get("dist_to_junction", 0)))[:4]

    return results, {"JLB": jlb, "JSB": jsb, "JSA": jsa, "JLA": jla}


def format_junction_report(variety, junctions, junction_pos, regions, genome_len):
    """Format junction analysis as text."""
    lines = []
    lines.append(f"\n{'='*60}")
    lines.append(f"{variety} (genome: {genome_len:,} bp)")
    lines.append(f"  LSC: 1..{regions['LSC'][1]:,} ({regions['LSC'][1]-regions['LSC'][0]:,} bp)")
    lines.append(f"  IRb: {regions['IRb'][0]+1:,}..{regions['IRb'][1]:,} ({regions['IRb'][1]-regions['IRb'][0]:,} bp)")
    lines.append(f"  SSC: {regions['SSC'][0]+1:,}..{regions['SSC'][1]:,} ({regions['SSC'][1]-regions['SSC'][0]:,} bp)")
    lines.append(f"  IRa: {regions['IRa'][0]+1:,}..{regions['IRa'][1]:,} ({regions['IRa'][1]-regions['IRa'][0]:,} bp)")

    for jname in ["JLB", "JSB", "JSA", "JLA"]:
        lines.append(f"\n  --- {jname} (position {junction_pos[jname]+1:,}) ---")
        for g in junctions[jname]:
            pos = g["position"]
            pseudo = " [pseudo]" if g.get("pseudo") else ""
            if pos == "spanning":
                extra_parts = []
                for k, v in g.items():
                    if k.startswith("bp_in_"):
                        extra_parts.append(f"{k}={v}")
                extra = ", ".join(extra_parts)
                lines.append(f"    {g['name']}{pseudo}: SPANNING junction ({extra})")
            else:
                lines.append(f"    {g['name']}{pseudo}: {pos}, {g['dist_to_junction']} bp from junction "
                           f"({g['start']+1:,}..{g['end']:,}, {g['length']} bp, strand={'+'if g['strand']==1 else '-'})")

    return "\n".join(lines)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("IR Boundary (IRscope) Analysis")
    print("=" * 40)

    all_lines = []
    csv_rows = []

    for var in VARIETIES:
        gb_path = os.path.join(GB_DIR, f"{var}.gb")
        record = SeqIO.read(gb_path, "genbank")
        genome_len = len(record.seq)

        regions = get_regions(record)
        genes = get_gene_features(record)
        junctions, junction_pos = find_junction_genes(genes, regions, genome_len)

        report = format_junction_report(var, junctions, junction_pos, regions, genome_len)
        all_lines.append(report)
        print(report)

        # CSV row for each junction
        for jname in ["JLB", "JSB", "JSA", "JLA"]:
            for g in junctions[jname]:
                row = {
                    "Variety": var,
                    "Junction": jname,
                    "Junction_pos": junction_pos[jname] + 1,
                    "Gene": g["name"],
                    "Pseudo": "yes" if g.get("pseudo") else "",
                    "Position": g["position"],
                    "Gene_start": g["start"] + 1,
                    "Gene_end": g["end"],
                    "Length": g["length"],
                    "Strand": "+" if g["strand"] == 1 else "-",
                    "Dist_to_junction": g.get("dist_to_junction", 0),
                }
                if g["position"] == "spanning":
                    for k, v in g.items():
                        if k.startswith("bp_in_"):
                            row[k] = v
                csv_rows.append(row)

    # Write text report
    report_path = os.path.join(OUT_DIR, "ir_boundary_report.txt")
    with open(report_path, "w") as f:
        f.write("IR Boundary (IRscope-style) Analysis\n")
        f.write("Rubus crataegifolius chloroplast genomes\n")
        f.write("\n".join(all_lines))
    print(f"\nText report: {report_path}")

    # Write CSV
    csv_path = os.path.join(OUT_DIR, "ir_boundary_genes.csv")
    fieldnames = ["Variety", "Junction", "Junction_pos", "Gene", "Pseudo",
                  "Position", "Gene_start", "Gene_end", "Length", "Strand",
                  "Dist_to_junction", "bp_in_LSC", "bp_in_IRb", "bp_in_SSC", "bp_in_IRa"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in csv_rows:
            writer.writerow(row)
    print(f"CSV: {csv_path}")

    # Write summary table (paper-ready)
    summary_path = os.path.join(OUT_DIR, "ir_boundary_summary.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Variety", "Genome_bp",
                         "LSC_bp", "IRb_bp", "SSC_bp", "IRa_bp",
                         "JLB_gene_LSC", "JLB_dist_LSC",
                         "JLB_gene_IRb", "JLB_dist_IRb",
                         "JSB_gene_IRb", "JSB_dist_IRb",
                         "JSB_gene_SSC", "JSB_dist_SSC",
                         "JSA_gene_SSC", "JSA_dist_SSC",
                         "JSA_gene_IRa", "JSA_dist_IRa",
                         "JLA_gene_IRa", "JLA_dist_IRa",
                         "JLA_gene_LSC", "JLA_dist_LSC"])

        for var in VARIETIES:
            gb_path = os.path.join(GB_DIR, f"{var}.gb")
            record = SeqIO.read(gb_path, "genbank")
            genome_len = len(record.seq)
            regions = get_regions(record)
            genes = get_gene_features(record)
            junctions, junction_pos = find_junction_genes(genes, regions, genome_len)

            lsc_len = regions["LSC"][1] - regions["LSC"][0]
            irb_len = regions["IRb"][1] - regions["IRb"][0]
            ssc_len = regions["SSC"][1] - regions["SSC"][0]
            ira_len = regions["IRa"][1] - regions["IRa"][0]

            def get_nearest(jname, side):
                """Get nearest gene on specified side of junction."""
                for g in junctions[jname]:
                    if g["position"] == side:
                        return g["name"] + (" [ψ]" if g.get("pseudo") else ""), g["dist_to_junction"]
                    elif g["position"] == "spanning":
                        key = f"bp_in_{side}"
                        if key in g:
                            return g["name"] + " (span)", g.get(key, 0)
                return "", ""

            jlb_lsc = get_nearest("JLB", "LSC")
            jlb_irb = get_nearest("JLB", "IRb")
            jsb_irb = get_nearest("JSB", "IRb")
            jsb_ssc = get_nearest("JSB", "SSC")
            jsa_ssc = get_nearest("JSA", "SSC")
            jsa_ira = get_nearest("JSA", "IRa")
            jla_ira = get_nearest("JLA", "IRa")
            jla_lsc = get_nearest("JLA", "LSC")

            writer.writerow([var, genome_len,
                           lsc_len, irb_len, ssc_len, ira_len,
                           jlb_lsc[0], jlb_lsc[1],
                           jlb_irb[0], jlb_irb[1],
                           jsb_irb[0], jsb_irb[1],
                           jsb_ssc[0], jsb_ssc[1],
                           jsa_ssc[0], jsa_ssc[1],
                           jsa_ira[0], jsa_ira[1],
                           jla_ira[0], jla_ira[1],
                           jla_lsc[0], jla_lsc[1]])

    print(f"Summary: {summary_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
