#!/usr/bin/env python3
"""
tbl_to_gb.py — Convert BankIt .tbl + .fsa → GenBank .gb for OGDRAW

Parses the NCBI feature table format and combines with FASTA sequence
to produce GenBank files suitable for OGDRAW genome map rendering.
"""

import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]
BANKIT_DIR = f"{BASE}/bankit"
OUT_DIR = f"{BASE}/bankit_gb"


def parse_tbl(tbl_path):
    """Parse NCBI feature table (.tbl) into a list of SeqFeature objects."""
    features = []
    current_feat = None
    location_parts = []

    with open(tbl_path) as f:
        for line in f:
            line = line.rstrip("\n")

            # Skip header line
            if line.startswith(">Feature"):
                continue

            # Check if this is a location/feature type line (starts with digit)
            match = re.match(r'^(\d+)\t(\d+)\t(\S+)', line)
            if match:
                # Save previous feature
                if current_feat is not None:
                    feat = _build_feature(current_feat, location_parts)
                    if feat:
                        features.append(feat)

                # Start new feature
                start, end, feat_type = match.groups()
                current_feat = {
                    "type": feat_type,
                    "qualifiers": {}
                }
                location_parts = [(int(start), int(end))]
                continue

            # Check if this is a continuation location line
            match = re.match(r'^(\d+)\t(\d+)\s*$', line)
            if match:
                start, end = match.groups()
                location_parts.append((int(start), int(end)))
                continue

            # Check if this is a qualifier line
            match = re.match(r'^\t\t\t(\S+)\t?(.*)', line)
            if match:
                key, val = match.groups()
                if current_feat is not None:
                    if key in current_feat["qualifiers"]:
                        current_feat["qualifiers"][key].append(val)
                    else:
                        current_feat["qualifiers"][key] = [val]
                continue

    # Save last feature
    if current_feat is not None:
        feat = _build_feature(current_feat, location_parts)
        if feat:
            features.append(feat)

    return features


def _build_feature(feat_dict, location_parts):
    """Build a BioPython SeqFeature from parsed data."""
    if not location_parts:
        return None

    locations = []
    for start, end in location_parts:
        if start > end:
            # Complement strand
            loc = FeatureLocation(end - 1, start, strand=-1)
        else:
            # Forward strand
            loc = FeatureLocation(start - 1, end, strand=1)
        locations.append(loc)

    if len(locations) == 1:
        location = locations[0]
    else:
        location = CompoundLocation(locations)

    # Clean up qualifiers
    quals = {}
    for key, vals in feat_dict["qualifiers"].items():
        if key == "pseudo":
            quals[key] = [""]
        elif key == "trans_splicing":
            quals[key] = [""]
        else:
            quals[key] = vals

    return SeqFeature(
        location=location,
        type=feat_dict["type"],
        qualifiers=quals
    )


def convert_variety(variety):
    """Convert one variety's .tbl + .fsa to .gb."""
    tbl_path = os.path.join(BANKIT_DIR, f"{variety}.tbl")
    fsa_path = os.path.join(BANKIT_DIR, f"{variety}.fsa")
    gb_path = os.path.join(OUT_DIR, f"{variety}.gb")

    # Read sequence
    fasta_rec = SeqIO.read(fsa_path, "fasta")
    seq = str(fasta_rec.seq)

    # Parse features
    features = parse_tbl(tbl_path)

    # Build GenBank record
    record = SeqRecord(
        Seq(seq),
        id=variety,
        name=variety,
        description=f"Rubus crataegifolius cultivar {variety} chloroplast, complete genome",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "data_file_division": "PLN",
            "organism": "Rubus crataegifolius",
            "source": "plastid:chloroplast",
        }
    )
    record.features = features

    # Write GenBank
    SeqIO.write(record, gb_path, "genbank")
    print(f"  {variety}: {len(features)} features → {gb_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Output: {OUT_DIR}\n")

    for variety in VARIETIES:
        try:
            convert_variety(variety)
        except Exception as e:
            print(f"  ERROR {variety}: {e}", file=sys.stderr)

    print(f"\nDone. Upload .gb files to OGDRAW (https://chlorobox.mpimp-golm.mpg.de/OGDraw.html)")


if __name__ == "__main__":
    main()
