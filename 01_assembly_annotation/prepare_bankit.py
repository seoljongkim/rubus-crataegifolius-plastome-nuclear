#!/usr/bin/env python3
"""
prepare_bankit.py — GeSeq .gb → BankIt submission files (.fsa + .tbl)

Processes 10 R. crataegifolius chloroplast genome GeSeq GenBank files:
1. Removes duplicate features (keeps best annotator per gene)
2. Fixes metadata (organism, cultivar, country)
3. Removes non-NCBI qualifiers (/info, /annotator)
4. Adds /transl_table=11 to CDS
5. Outputs .fsa (FASTA) + .tbl (Feature Table) per variety
"""

import os
import re
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# ── Config ──────────────────────────────────────────────────────────────────

VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]

GB_DIR = f"{BASE}/job-results-20263552725"
GB_TEMPLATE = "GeSeqJob-20260304-103051_{variety}_GenBank.gb"

FASTA_DIR = f"{BASE}/standardized"
FASTA_TEMPLATE = "{variety}_std.fasta"

OUT_DIR = f"{BASE}/bankit"

# Annotator priority (higher index = higher priority)
ANNOTATOR_PRIORITY = {
    "HMMER": 1,
    "blatX": 2,
    "blatN": 3,
    "ARAGORN": 4,
    "tRNAscan-SE_v2": 5,
    "Chloe": 6,
    "merged": 7,
}

# Feature types that represent a gene's functional annotation
GENE_FEATURE_TYPES = {"gene", "CDS", "tRNA", "rRNA", "exon", "intron"}

# Gene name corrections: GeSeq pipe-separated names → standard names
# Based on NC_039704.1 reference (blatN hits) + Rubus literature convention
GENE_NAME_FIXES = {
    "trnG-GCC | trnG-UCC": "trnG-UCC",
    "trnG-UCC | trnG-GCC": "trnG-UCC",
    "trnI-GAU | trnE-UUC | trnl-GAU": "trnI-GAU",
    "trnI-GAU | trnl-GAU": "trnI-GAU",
    "trnL-UAA | trnL": "trnL-UAA",
    "trnM-CAU | trnI-CAU": "trnI-CAU",
    "trnM-CAU | trnfM-CAU": "trnfM-CAU",
    "trnV-UAC | trnI-AAU": "trnV-UAC",
    # GeSeq Chloe new names → traditional names (Rubus literature standard)
    "pafI": "ycf3",
    "pafII": "ycf4",
    "clpP1": "clpP",
    "pbf1": "psbN",
}

# Standalone gene names to exclude (spurious ARAGORN-only annotations)
EXCLUDE_GENES = {"trnG-GCC"}

# tRNA product mapping for corrected gene names
TRNA_PRODUCT = {
    "trnG-UCC": "tRNA-Gly",
    "trnI-GAU": "tRNA-Ile",
    "trnL-UAA": "tRNA-Leu",
    "trnI-CAU": "tRNA-Ile",
    "trnfM-CAU": "tRNA-fMet",
    "trnV-UAC": "tRNA-Val",
}

# Qualifiers to remove (not NCBI-compatible)
REMOVE_QUALIFIERS = {"info", "annotator"}

# Qualifiers to keep
KEEP_QUALIFIERS = {
    "gene", "product", "translation", "codon_start", "trans_splicing",
    "note", "transl_table", "rpt_type", "organelle", "organism",
    "mol_type", "cultivar", "country", "pseudo",
}


# ── Helpers ─────────────────────────────────────────────────────────────────

def parse_annotator(feature):
    """Extract annotator string from /annotator qualifier."""
    annots = feature.qualifiers.get("annotator", [""])
    return annots[0] if annots else ""


def annotator_score(annotator_str):
    """Score an annotator string. Higher = better priority."""
    if not annotator_str:
        return 0
    # Check for "merged" keyword
    if "merged" in annotator_str:
        # Check if it's a plain merged or specific tool merged
        return ANNOTATOR_PRIORITY["merged"]
    # Find best individual annotator
    best = 0
    for tool, score in ANNOTATOR_PRIORITY.items():
        if tool in annotator_str:
            best = max(best, score)
    return best


def feature_span(feature):
    """Total span of a feature in bp."""
    return sum(len(part) for part in feature.location.parts)


def gene_name(feature):
    """Extract gene name from /gene qualifier."""
    names = feature.qualifiers.get("gene", [])
    return names[0] if names else None


def fix_gene_name(name):
    """Fix pipe-separated gene names from GeSeq."""
    if name and name in GENE_NAME_FIXES:
        return GENE_NAME_FIXES[name]
    return name


def location_key(feature):
    """Create a hashable key from feature location for grouping."""
    parts = []
    for part in feature.location.parts:
        parts.append((int(part.start), int(part.end), part.strand))
    return tuple(sorted(parts))


def deduplicate_features(record):
    """
    Remove duplicate features, keeping the best annotator per gene per type.
    Returns a list of cleaned features.
    """
    # Separate structural features from gene features
    structural = []  # source, misc_feature, repeat_region
    gene_features = defaultdict(list)  # gene_name -> list of features

    for feat in record.features:
        if feat.type in ("source", "misc_feature", "repeat_region"):
            structural.append(feat)
            continue
        gname = gene_name(feat)
        if gname:
            # Fix gene name first (pipe-separated names)
            fixed = fix_gene_name(gname)
            # Skip excluded genes (e.g., standalone trnG-GCC)
            if fixed in EXCLUDE_GENES:
                continue
            gene_features[gname].append(feat)
        else:
            # Features without gene name — keep as-is (rare)
            structural.append(feat)

    # Process each gene group
    result_features = list(structural)

    for gname, feats in gene_features.items():
        selected = select_best_features(gname, feats)
        result_features.extend(selected)

    # Add ycf1 pseudogene at IRb/SSC boundary
    ycf1_pseudo = _make_ycf1_pseudogene(result_features)
    if ycf1_pseudo:
        result_features.append(ycf1_pseudo)

    # Ensure IR-duplicated genes have annotations in both IRb and IRa
    result_features = _ensure_ir_duplicates(result_features)

    # Sort by start position, then by feature type priority
    type_order = {"source": 0, "misc_feature": 1, "repeat_region": 2,
                  "gene": 3, "mRNA": 4, "CDS": 5, "tRNA": 6, "rRNA": 7,
                  "exon": 8, "intron": 9}
    result_features.sort(key=lambda f: (
        min(int(p.start) for p in f.location.parts),
        type_order.get(f.type, 99)
    ))

    return result_features


def _make_ycf1_pseudogene(features):
    """
    Create a ycf1 pseudogene gene feature at the IRb/SSC boundary.

    The full ycf1 spans SSC→IRa on the complement strand. The truncated
    copy in IRb is a pseudogene (no functional CDS). We calculate its
    position by mirroring the IRa portion of ycf1 into IRb.
    """
    # Find IRb end position (from repeat_region for IRB)
    irb_end = None
    ira_start = None
    for f in features:
        if f.type == "repeat_region":
            note = f.qualifiers.get("note", [""])[0]
            if "IRB" in note or "inverted repeat B" in note.lower():
                irb_end = int(f.location.end)
            elif "IRA" in note or "inverted repeat A" in note.lower():
                ira_start = int(f.location.start)

    # Find full ycf1 gene position
    ycf1_end = None
    for f in features:
        if f.type == "gene" and gene_name(f) == "ycf1":
            ycf1_end = max(int(p.end) for p in f.location.parts)

    if irb_end is None or ira_start is None or ycf1_end is None:
        return None

    # Calculate overlap length of ycf1 in IRa
    overlap_in_ira = ycf1_end - ira_start
    if overlap_in_ira <= 0:
        return None

    # Pseudogene in IRb: complement strand, mirroring the IRa portion
    pseudo_start = irb_end - overlap_in_ira  # 0-based start
    pseudo_end = irb_end                      # 0-based exclusive end

    pseudo = SeqFeature(
        FeatureLocation(pseudo_start, pseudo_end, strand=-1),
        type="gene",
        qualifiers={
            "gene": ["ycf1"],
            "pseudo": [""],
            "note": ["truncated copy of ycf1 in IRb"],
        }
    )
    return pseudo


# IR-duplicated genes (present in both IRb and IRa)
# Excluding rps12 (trans-spliced, handled separately)
IR_DUPLICATED_GENES = {
    "rpl2", "rpl23", "trnI-CAU", "ycf2", "trnL-CAA", "ndhB",
    "rps7", "trnV-GAC", "rrn16", "trnI-GAU", "trnA-UGC",
    "rrn23", "rrn4.5", "rrn5", "trnR-ACG", "trnN-GUU",
}


def _find_ir_boundaries(features):
    """Find IRb and IRa boundaries from repeat_region features. Returns 0-based."""
    irb_start = irb_end = ira_start = ira_end = None
    for f in features:
        if f.type == "repeat_region":
            note = f.qualifiers.get("note", [""])[0]
            if "IRB" in note or "inverted repeat B" in note.lower():
                irb_start = int(f.location.start)
                irb_end = int(f.location.end)
            elif "IRA" in note or "inverted repeat A" in note.lower():
                ira_start = int(f.location.start)
                ira_end = int(f.location.end)
    return irb_start, irb_end, ira_start, ira_end


def _feature_ir_region(feat, irb_start, irb_end, ira_start, ira_end):
    """Determine which IR region a feature is primarily in."""
    # Use the midpoint of the feature (excluding LSC parts for trans-spliced)
    for part in feat.location.parts:
        mid = (int(part.start) + int(part.end)) / 2
        if irb_start <= mid <= irb_end:
            return "IRb"
        if ira_start <= mid <= ira_end:
            return "IRa"
    return "other"


def _mirror_feature(feat, M):
    """
    Create a mirror copy of a feature in the opposite IR.
    M = irb_start + ira_end (0-based mirror constant).
    For position X in one IR, mirror = M - X.
    Strand flips: +1 → -1, -1 → +1.
    """
    new_parts = []
    for part in feat.location.parts:
        new_start = M - int(part.end)
        new_end = M - int(part.start)
        new_strand = -part.strand if part.strand else part.strand
        new_parts.append(FeatureLocation(new_start, new_end, strand=new_strand))

    if len(new_parts) == 1:
        new_location = new_parts[0]
    else:
        new_location = CompoundLocation(new_parts)

    return SeqFeature(
        location=new_location,
        type=feat.type,
        qualifiers=dict(feat.qualifiers),
    )


def _ensure_ir_duplicates(features):
    """
    Ensure IR-duplicated genes have complete annotations in both IRb and IRa.

    After GeSeq dedup, IR genes typically have:
    - gene feature in IRa only
    - CDS/tRNA/rRNA feature in IRb only
    This function adds the missing mirror copies, and also:
    - Adds a second rps12 CDS with exons 2-3 mirrored to the other IR
    - Adds trnL-CAG in IRa
    """
    irb_start, irb_end, ira_start, ira_end = _find_ir_boundaries(features)
    if irb_start is None:
        return features

    M = irb_start + ira_end  # mirror constant (0-based)

    # Collect IR gene features and non-IR features separately
    ir_gene_feats = defaultdict(list)
    other_features = []

    for f in features:
        gn = gene_name(f)
        if gn in IR_DUPLICATED_GENES:
            ir_gene_feats[gn].append(f)
        else:
            other_features.append(f)

    result = list(other_features)

    for gn, feats in ir_gene_feats.items():
        # Separate features by IR region
        irb_feats = []
        ira_feats = []
        for f in feats:
            region = _feature_ir_region(f, irb_start, irb_end, ira_start, ira_end)
            if region == "IRb":
                irb_feats.append(f)
            elif region == "IRa":
                ira_feats.append(f)

        # Keep only essential feature types: gene, CDS/tRNA/rRNA
        # (skip redundant exon/intron from secondary annotators)
        irb_essential = [f for f in irb_feats
                         if f.type in ("gene", "CDS", "tRNA", "rRNA")]
        ira_essential = [f for f in ira_feats
                         if f.type in ("gene", "CDS", "tRNA", "rRNA")]

        irb_types = {f.type for f in irb_essential}
        ira_types = {f.type for f in ira_essential}

        # Start with all essential features from both IRs
        result.extend(irb_essential)
        result.extend(ira_essential)

        # Add missing mirrors
        for f in irb_essential:
            if f.type not in ira_types:
                result.append(_mirror_feature(f, M))
        for f in ira_essential:
            if f.type not in irb_types:
                result.append(_mirror_feature(f, M))

    # Handle rps12 second copy (trans-spliced)
    rps12_mirror = _make_rps12_ir_copy(features, irb_start, irb_end,
                                        ira_start, ira_end, M)
    result.extend(rps12_mirror)

    # Add trnL-CAG in IRa
    trnl_cag = _make_trnl_cag(features, ira_start, ira_end, M)
    result.extend(trnl_cag)

    return result


def _make_rps12_ir_copy(features, irb_start, irb_end, ira_start, ira_end, M):
    """
    Create the second rps12 copy with exons 2-3 in the other IR.

    rps12 is trans-spliced: exon 1 in LSC, exons 2-3 in one IR.
    We find the existing CDS, identify the IR exons, mirror them, and
    create a second CDS + gene with mirrored exons + same LSC exon 1.
    """
    result = []

    # Find existing rps12 CDS
    rps12_cds = None
    for f in features:
        if f.type == "CDS" and gene_name(f) == "rps12":
            rps12_cds = f
            break

    if rps12_cds is None:
        return result

    # Separate parts into LSC (exon 1) and IR (exons 2-3)
    lsc_parts = []
    ir_parts = []
    for part in rps12_cds.location.parts:
        mid = (int(part.start) + int(part.end)) / 2
        if irb_start <= mid <= irb_end or ira_start <= mid <= ira_end:
            ir_parts.append(part)
        else:
            lsc_parts.append(part)

    if not lsc_parts or not ir_parts:
        return result

    # Mirror IR parts to the other IR
    mirrored_ir_parts = []
    for part in ir_parts:
        new_start = M - int(part.end)
        new_end = M - int(part.start)
        new_strand = -part.strand if part.strand else part.strand
        mirrored_ir_parts.append(FeatureLocation(new_start, new_end, strand=new_strand))

    # Build the new CDS: LSC exon(s) + mirrored IR exons
    # Order: for the mirrored copy, the exons should be ordered properly
    all_parts = lsc_parts + mirrored_ir_parts

    new_cds = SeqFeature(
        location=CompoundLocation(all_parts),
        type="CDS",
        qualifiers=dict(rps12_cds.qualifiers),
    )
    result.append(new_cds)

    # Build gene feature spanning all parts
    all_positions = []
    for part in all_parts:
        all_positions.extend([int(part.start), int(part.end)])
    gene_start = min(all_positions)
    gene_end = max(all_positions)

    # Gene feature needs to encompass all parts; use CompoundLocation if disjoint
    new_gene = SeqFeature(
        location=CompoundLocation(all_parts),
        type="gene",
        qualifiers={
            "gene": ["rps12"],
            "trans_splicing": [""],
        },
    )
    result.append(new_gene)

    return result


def _make_trnl_cag(features, ira_start, ira_end, M):
    """
    Add trnL-CAG in IRa.

    trnL-CAG is located between trnR-ACG and trnN-GUU in IRa.
    We compute its position relative to trnR-ACG.
    In DAHO (CP_ref): trnR-ACG IRa = complement(131831..131904),
                       trnL-CAG IRa = 131309..131370
    Offset from trnR-ACG start (131830, 0-based): 131308 - 131830 = -522
    ... actually let's compute it from the mirror constant.
    In DAHO: trnR-ACG in IRb = 109135..109208 (+ strand, 0-based: 109134..109208)
    trnL-CAG should be at IRa position ~131309..131370 (1-based) = 131308..131370 (0-based)
    trnR-ACG in IRa: complement(131831..131904) (1-based) = 131830..131904 (0-based), -strand
    Distance from trnR-ACG IRa end to trnL-CAG start: 131830 - 131370 = 460 (0-based)
    """
    # Find trnR-ACG in IRa to compute trnL-CAG position
    trnr_ira = None
    for f in features:
        if gene_name(f) == "trnR-ACG":
            region = _feature_ir_region(f, M - ira_end, M - ira_start,
                                        ira_start, ira_end)
            # Check if this feature is in IRa
            for part in f.location.parts:
                mid = (int(part.start) + int(part.end)) / 2
                if ira_start <= mid <= ira_end:
                    trnr_ira = f
                    break
            if trnr_ira:
                break

    if trnr_ira is None:
        return []

    # trnL-CAG is 62 bp, located before trnR-ACG in IRa
    # In DAHO: trnR-ACG IRa = complement(131831..131904) → 0-based [131830, 131904)
    # trnL-CAG IRa = 131309..131370 → 0-based [131308, 131370), + strand
    # Gap between trnL-CAG end and trnR-ACG start: 131830 - 131370 = 460

    trnr_start = min(int(p.start) for p in trnr_ira.location.parts)
    trnl_end = trnr_start - 460  # 0-based end (exclusive)
    trnl_start = trnl_end - 62   # 0-based start

    # Verify it's within IRa
    if trnl_start < ira_start or trnl_end > ira_end:
        return []

    result = []
    result.append(SeqFeature(
        FeatureLocation(trnl_start, trnl_end, strand=1),
        type="gene",
        qualifiers={"gene": ["trnL-CAG"]},
    ))
    result.append(SeqFeature(
        FeatureLocation(trnl_start, trnl_end, strand=1),
        type="tRNA",
        qualifiers={"gene": ["trnL-CAG"], "product": ["tRNA-Leu"]},
    ))
    return result


def select_best_features(gname, feats):
    """
    For a given gene, select the best annotator's features.
    Returns: one gene feature + one functional feature (CDS/tRNA/rRNA)
             + associated exon/intron features.
    """
    # Group by feature type
    by_type = defaultdict(list)
    for f in feats:
        by_type[f.type].append(f)

    selected = []

    # 1. Select best functional feature (CDS, tRNA, or rRNA)
    func_type = None
    for ft in ("CDS", "tRNA", "rRNA"):
        if ft in by_type:
            func_type = ft
            break

    best_func = None
    best_annotator_str = ""
    if func_type:
        candidates = by_type[func_type]
        # Special case: atpF — prefer "HMMER, blatX; merged" over Chloe
        if gname == "atpF":
            best_func = _select_atpf_best(candidates)
        else:
            best_func = max(candidates, key=lambda f: (
                annotator_score(parse_annotator(f)),
                feature_span(f)
            ))
        if best_func:
            best_annotator_str = parse_annotator(best_func)
            selected.append(best_func)

    # 2. Select best gene feature (widest span from best annotator)
    if "gene" in by_type:
        gene_feats = by_type["gene"]
        # Prefer gene feature with matching annotator, or widest span
        best_gene = None
        for gf in gene_feats:
            ann = parse_annotator(gf)
            if best_gene is None:
                best_gene = gf
            elif _annotators_compatible(ann, best_annotator_str):
                if best_gene is None or not _annotators_compatible(
                        parse_annotator(best_gene), best_annotator_str):
                    best_gene = gf
                elif feature_span(gf) >= feature_span(best_gene):
                    best_gene = gf
            elif not _annotators_compatible(
                    parse_annotator(best_gene), best_annotator_str):
                if feature_span(gf) > feature_span(best_gene):
                    best_gene = gf
        if best_gene:
            selected.append(best_gene)

    # 3. Select exon/intron features matching the best annotator
    for ft in ("exon", "intron"):
        if ft in by_type:
            for f in by_type[ft]:
                ann = parse_annotator(f)
                if _annotators_compatible(ann, best_annotator_str):
                    selected.append(f)
                    break  # one per exon/intron number? Actually keep all matching

    # Actually, for exons/introns we want all that match the annotator
    # Re-do exon/intron selection properly
    selected = [s for s in selected if s.type not in ("exon", "intron")]
    for ft in ("exon", "intron"):
        if ft in by_type:
            matching = [f for f in by_type[ft]
                        if _annotators_compatible(parse_annotator(f), best_annotator_str)]
            if matching:
                selected.extend(matching)
            # If no matching annotator for exon/intron, try the best-scored ones
            elif by_type[ft]:
                best_ann = max(by_type[ft], key=lambda f: annotator_score(parse_annotator(f)))
                target_ann = parse_annotator(best_ann)
                selected.extend(f for f in by_type[ft]
                                if parse_annotator(f) == target_ann)

    return selected


def _select_atpf_best(candidates):
    """For atpF, prefer 555 bp CDS (HMMER/blatX merged), not Chloe 145 bp."""
    # Prefer the longest CDS that is NOT purely Chloe
    non_chloe = [c for c in candidates if "Chloe" not in parse_annotator(c)
                 or "merged" in parse_annotator(c)]
    if non_chloe:
        # Among non-Chloe, pick longest
        return max(non_chloe, key=lambda f: feature_span(f))
    return max(candidates, key=lambda f: feature_span(f))


def _annotators_compatible(ann1, ann2):
    """Check if two annotator strings share any tool."""
    if not ann1 or not ann2:
        return False
    tools1 = set(re.split(r'[,;\s]+', ann1))
    tools2 = set(re.split(r'[,;\s]+', ann2))
    return bool(tools1 & tools2)


def clean_qualifiers(feature, variety):
    """Remove non-NCBI qualifiers, fix gene names, add /transl_table=11 for CDS."""
    new_quals = {}
    for key, val in feature.qualifiers.items():
        if key in REMOVE_QUALIFIERS:
            continue
        new_quals[key] = val

    # Fix pipe-separated gene names
    if "gene" in new_quals:
        original = new_quals["gene"][0]
        fixed = fix_gene_name(original)
        if fixed != original:
            new_quals["gene"] = [fixed]
            # Also fix tRNA product if needed
            if "product" in new_quals and fixed in TRNA_PRODUCT:
                new_quals["product"] = [TRNA_PRODUCT[fixed]]

    # Add transl_table for CDS
    if feature.type == "CDS":
        new_quals["transl_table"] = ["11"]

    # Fix source metadata
    if feature.type == "source":
        new_quals = {
            "organism": ["Rubus crataegifolius"],
            "organelle": ["plastid:chloroplast"],
            "mol_type": ["genomic DNA"],
            "cultivar": [variety],
            "country": ["South Korea"],
        }

    feature.qualifiers = new_quals
    return feature


# ── .tbl output ─────────────────────────────────────────────────────────────

def write_tbl(features, seq_id, genome_size, out_path):
    """Write NCBI feature table format."""
    with open(out_path, "w") as f:
        f.write(f">Feature {seq_id}\n")
        for feat in features:
            _write_feature_tbl(f, feat, genome_size)


def _write_feature_tbl(f, feat, genome_size):
    """Write a single feature in .tbl format."""
    parts = list(feat.location.parts)

    if len(parts) == 0:
        return

    # Write location lines
    for i, part in enumerate(parts):
        start = int(part.start) + 1  # Convert to 1-based
        end = int(part.end)

        if part.strand == -1:
            # Complement: swap start/end
            start, end = end, start

        if i == 0:
            f.write(f"{start}\t{end}\t{feat.type}\n")
        else:
            f.write(f"{start}\t{end}\n")

    # Write qualifiers
    for key, values in feat.qualifiers.items():
        if key in REMOVE_QUALIFIERS:
            continue
        for val in values:
            if val is True or val == "":
                # Boolean qualifier (e.g., /trans_splicing)
                f.write(f"\t\t\t{key}\n")
            else:
                f.write(f"\t\t\t{key}\t{val}\n")


# ── .fsa output ─────────────────────────────────────────────────────────────

def write_fsa(seq_str, variety, genome_size, out_path):
    """Write BankIt FASTA with proper header."""
    header = (
        f">{variety} "
        f"[organism=Rubus crataegifolius] "
        f"[organelle=plastid:chloroplast] "
        f"[moltype=genomic DNA] "
        f"[topology=circular] "
        f"[cultivar={variety}] "
        f"[country=South Korea]"
    )
    with open(out_path, "w") as f:
        f.write(header + "\n")
        # Write sequence in 70-char lines
        seq = str(seq_str).upper()
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")


# ── Main ────────────────────────────────────────────────────────────────────

def process_variety(variety):
    """Process one variety: read .gb, deduplicate, write .tbl + .fsa."""
    gb_path = os.path.join(GB_DIR, GB_TEMPLATE.format(variety=variety))
    fasta_path = os.path.join(FASTA_DIR, FASTA_TEMPLATE.format(variety=variety))

    if not os.path.exists(gb_path):
        print(f"  ERROR: {gb_path} not found", file=sys.stderr)
        return False
    if not os.path.exists(fasta_path):
        print(f"  ERROR: {fasta_path} not found", file=sys.stderr)
        return False

    # Parse GenBank
    record = SeqIO.read(gb_path, "genbank")

    # Read standardized FASTA sequence
    fasta_record = SeqIO.read(fasta_path, "fasta")
    seq_str = str(fasta_record.seq)
    genome_size = len(seq_str)

    # Verify sequences match
    gb_seq = str(record.seq).upper()
    if gb_seq != seq_str.upper():
        print(f"  WARNING: {variety} sequences differ! "
              f"GB={len(gb_seq)} bp, FASTA={genome_size} bp", file=sys.stderr)

    print(f"  {variety}: {genome_size} bp, {len(record.features)} raw features")

    # Deduplicate features
    features = deduplicate_features(record)
    print(f"  → {len(features)} features after dedup")

    # Clean qualifiers
    for feat in features:
        clean_qualifiers(feat, variety)

    # Count gene features for reporting
    gene_count = sum(1 for f in features if f.type == "gene")
    cds_count = sum(1 for f in features if f.type == "CDS")
    trna_count = sum(1 for f in features if f.type == "tRNA")
    rrna_count = sum(1 for f in features if f.type == "rRNA")
    print(f"  → genes={gene_count}, CDS={cds_count}, tRNA={trna_count}, rRNA={rrna_count}")

    # Check atpF CDS size
    for feat in features:
        if feat.type == "CDS" and gene_name(feat) == "atpF":
            span = feature_span(feat)
            print(f"  → atpF CDS: {span} bp {'✓' if span == 555 else '✗ UNEXPECTED'}")
            break

    # Check rps18 CDS size for GB048
    if variety == "GB048":
        for feat in features:
            if feat.type == "CDS" and gene_name(feat) == "rps18":
                span = feature_span(feat)
                print(f"  → rps18 CDS: {span} bp {'✓' if span == 312 else '✗ UNEXPECTED'}")
                break

    # Check rps12 trans_splicing
    rps12_ts = any(
        feat.type == "CDS" and gene_name(feat) == "rps12"
        and "trans_splicing" in feat.qualifiers
        for feat in features
    )
    print(f"  → rps12 trans_splicing: {'✓' if rps12_ts else '✗ MISSING'}")

    # Write output files
    tbl_path = os.path.join(OUT_DIR, f"{variety}.tbl")
    fsa_path = os.path.join(OUT_DIR, f"{variety}.fsa")

    write_tbl(features, variety, genome_size, tbl_path)
    write_fsa(seq_str, variety, genome_size, fsa_path)

    print(f"  → {tbl_path}")
    print(f"  → {fsa_path}")
    return True


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Output directory: {OUT_DIR}\n")

    success = 0
    for variety in VARIETIES:
        print(f"Processing {variety}...")
        if process_variety(variety):
            success += 1
        print()

    print(f"Done: {success}/{len(VARIETIES)} varieties processed.")

    # Summary check
    print("\n=== Output files ===")
    for variety in VARIETIES:
        tbl = os.path.join(OUT_DIR, f"{variety}.tbl")
        fsa = os.path.join(OUT_DIR, f"{variety}.fsa")
        tbl_ok = "✓" if os.path.exists(tbl) else "✗"
        fsa_ok = "✓" if os.path.exists(fsa) else "✗"
        print(f"  {variety}: .tbl {tbl_ok}  .fsa {fsa_ok}")


if __name__ == "__main__":
    main()
