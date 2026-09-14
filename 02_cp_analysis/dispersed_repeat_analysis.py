#!/usr/bin/env python3
"""
Dispersed repeat analysis for chloroplast genomes using BLASTn self-alignment.
Equivalent to REPuter with parameters: -l 30 -h 3 -f -p -r -c

Approach:
- BLASTn each genome against itself
- Filter for alignments >= 30 bp
- Categorize repeats as Forward (F), Palindromic (P), Reverse (R), Complement (C)
- Deduplicate: each repeat pair counted once
- Flag large IR-derived palindromic repeats separately

Author: Claude Code analysis pipeline
Date: 2026-03-09
"""

import os
import sys
import subprocess
import tempfile
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# ── Parameters ──────────────────────────────────────────────────────────────
MIN_REPEAT_LEN = 30
MAX_MISMATCHES = 3  # Hamming distance
VARIETIES = ["DAHE", "DAHO", "DB", "DS", "GB048", "GB066", "GDC", "GDM", "GDW", "HGB"]
INPUT_DIR = f"{BASE}/standardized"
OUTPUT_DIR = f"{BASE}/reputer_output"
BLAST_CMD = "blastn"

# Size bins for distribution
SIZE_BINS = [(30, 44), (45, 59), (60, 89), (90, float('inf'))]
SIZE_BIN_LABELS = ["30-44", "45-59", "60-89", ">=90"]

# IR size threshold: repeats larger than this in palindromic category are IR-derived
IR_THRESHOLD = 1000  # bp


def run_blast_self(fasta_path, variety, tmpdir):
    """Run BLASTn of genome against itself, return parsed hits."""
    # Create blast database
    db_path = os.path.join(tmpdir, f"{variety}_db")
    subprocess.run(
        [BLAST_CMD.replace("blastn", "makeblastdb"), "-in", fasta_path,
         "-dbtype", "nucl", "-out", db_path],
        check=True, capture_output=True
    )

    # We need to search in 4 modes:
    # Forward (F): plus/plus matches (same strand)
    # Palindromic (P): plus/minus matches (reverse complement)
    # For Reverse (R) and Complement (C), we need additional queries

    # Read the sequence
    record = SeqIO.read(fasta_path, "fasta")
    seq = str(record.seq)
    seq_len = len(seq)

    all_repeats = []

    # ── Run blastn: finds Forward (plus/plus) and Palindromic (plus/minus) ──
    outfile = os.path.join(tmpdir, f"{variety}_blast.txt")
    # Use outfmt 6 with custom fields
    # qseqid sseqid qstart qend sstart send qlen slen length nident mismatch pident qseq sseq sstrand
    blast_cmd = [
        BLAST_CMD, "-query", fasta_path, "-db", db_path,
        "-outfmt", "6 qseqid sseqid qstart qend sstart send length nident mismatch pident sstrand",
        "-evalue", "1000",
        "-word_size", "7",
        "-dust", "no",
        "-soft_masking", "false",
        "-max_target_seqs", "50000",
        "-max_hsps", "50000",
        "-num_threads", "4",
        "-strand", "both",
        "-out", outfile,
        "-ungapped"  # ungapped for Hamming distance matching
    ]

    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Try without -ungapped if it fails
        blast_cmd.remove("-ungapped")
        subprocess.run(blast_cmd, check=True, capture_output=True, text=True)

    # Parse BLAST output
    with open(outfile) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            qstart = int(fields[2])
            qend = int(fields[3])
            sstart = int(fields[4])
            send = int(fields[5])
            length = int(fields[6])
            nident = int(fields[7])
            mismatch = int(fields[8])
            pident = float(fields[9])
            strand = fields[10]

            # Filter by minimum length
            if length < MIN_REPEAT_LEN:
                continue

            # Filter by max mismatches (Hamming distance <= 3)
            if mismatch > MAX_MISMATCHES:
                continue

            # Determine repeat type based on strand
            if strand == "plus":
                repeat_type = "F"  # Forward
            else:
                repeat_type = "P"  # Palindromic (reverse complement)

            all_repeats.append({
                'type': repeat_type,
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'length': length,
                'nident': nident,
                'mismatch': mismatch,
                'pident': pident,
                'seq_len': seq_len
            })

    # ── For Reverse (R) and Complement (C) repeats ──
    # Reverse: query matches the reverse (not complement) of subject
    # Complement: query matches the complement (not reverse) of subject
    # These are rare in biological sequences but REPuter reports them.
    #
    # To find these, we create modified versions of the sequence:
    # - Reverse of sequence (reversed, not complemented)
    # - Complement of sequence (complemented, not reversed)

    # Create reverse sequence (not complemented)
    rev_seq = seq[::-1]
    rev_fasta = os.path.join(tmpdir, f"{variety}_rev.fasta")
    with open(rev_fasta, 'w') as f:
        f.write(f">{variety}_reverse\n{rev_seq}\n")

    # Create complement sequence (not reversed)
    comp_table = str.maketrans("ATCGatcg", "TAGCtagc")
    comp_seq = seq.translate(comp_table)
    comp_fasta = os.path.join(tmpdir, f"{variety}_comp.fasta")
    with open(comp_fasta, 'w') as f:
        f.write(f">{variety}_complement\n{comp_seq}\n")

    # Search for Reverse repeats: original query vs reversed-sequence database
    rev_db = os.path.join(tmpdir, f"{variety}_rev_db")
    subprocess.run(
        [BLAST_CMD.replace("blastn", "makeblastdb"), "-in", rev_fasta,
         "-dbtype", "nucl", "-out", rev_db],
        check=True, capture_output=True
    )

    rev_outfile = os.path.join(tmpdir, f"{variety}_rev_blast.txt")
    rev_blast_cmd = [
        BLAST_CMD, "-query", fasta_path, "-db", rev_db,
        "-outfmt", "6 qseqid sseqid qstart qend sstart send length nident mismatch pident sstrand",
        "-evalue", "1000",
        "-word_size", "7",
        "-dust", "no",
        "-soft_masking", "false",
        "-max_target_seqs", "50000",
        "-max_hsps", "50000",
        "-num_threads", "4",
        "-strand", "plus",  # only plus/plus since the db is already reversed
        "-out", rev_outfile,
        "-ungapped"
    ]

    result = subprocess.run(rev_blast_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        rev_blast_cmd.remove("-ungapped")
        subprocess.run(rev_blast_cmd, check=True, capture_output=True, text=True)

    with open(rev_outfile) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            qstart = int(fields[2])
            qend = int(fields[3])
            sstart = int(fields[4])
            send = int(fields[5])
            length = int(fields[6])
            nident = int(fields[7])
            mismatch = int(fields[8])
            pident = float(fields[9])

            if length < MIN_REPEAT_LEN:
                continue
            if mismatch > MAX_MISMATCHES:
                continue

            # Convert reversed coordinates back to original
            # In the reversed sequence, position p corresponds to original position (seq_len - p + 1)
            orig_sstart = seq_len - sstart + 1
            orig_send = seq_len - send + 1

            all_repeats.append({
                'type': 'R',
                'qstart': qstart,
                'qend': qend,
                'sstart': min(orig_sstart, orig_send),
                'send': max(orig_sstart, orig_send),
                'length': length,
                'nident': nident,
                'mismatch': mismatch,
                'pident': pident,
                'seq_len': seq_len
            })

    # Search for Complement repeats: original query vs complemented-sequence database
    comp_db = os.path.join(tmpdir, f"{variety}_comp_db")
    subprocess.run(
        [BLAST_CMD.replace("blastn", "makeblastdb"), "-in", comp_fasta,
         "-dbtype", "nucl", "-out", comp_db],
        check=True, capture_output=True
    )

    comp_outfile = os.path.join(tmpdir, f"{variety}_comp_blast.txt")
    comp_blast_cmd = [
        BLAST_CMD, "-query", fasta_path, "-db", comp_db,
        "-outfmt", "6 qseqid sseqid qstart qend sstart send length nident mismatch pident sstrand",
        "-evalue", "1000",
        "-word_size", "7",
        "-dust", "no",
        "-soft_masking", "false",
        "-max_target_seqs", "50000",
        "-max_hsps", "50000",
        "-num_threads", "4",
        "-strand", "plus",  # only plus/plus since the db is already complemented
        "-out", comp_outfile,
        "-ungapped"
    ]

    result = subprocess.run(comp_blast_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        comp_blast_cmd.remove("-ungapped")
        subprocess.run(comp_blast_cmd, check=True, capture_output=True, text=True)

    with open(comp_outfile) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            qstart = int(fields[2])
            qend = int(fields[3])
            sstart = int(fields[4])
            send = int(fields[5])
            length = int(fields[6])
            nident = int(fields[7])
            mismatch = int(fields[8])
            pident = float(fields[9])

            if length < MIN_REPEAT_LEN:
                continue
            if mismatch > MAX_MISMATCHES:
                continue

            all_repeats.append({
                'type': 'C',
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'length': length,
                'nident': nident,
                'mismatch': mismatch,
                'pident': pident,
                'seq_len': seq_len
            })

    return all_repeats, seq_len


def remove_self_hits(repeats):
    """Remove self-hits (where query and subject positions are identical on same strand)."""
    filtered = []
    for r in repeats:
        if r['type'] == 'F':
            # For forward repeats, self-hit is when positions are identical
            if r['qstart'] == r['sstart'] and r['qend'] == r['send']:
                continue
        elif r['type'] == 'P':
            # For palindromic, self-hit at same position on opposite strand
            # This can happen for palindromic sequences at the exact same location
            if r['qstart'] == r['sstart'] and r['qend'] == r['send']:
                continue
        filtered.append(r)
    return filtered


def deduplicate_repeats(repeats):
    """
    Deduplicate repeat pairs. For Forward repeats, (A,B) and (B,A) are the same pair.
    For Palindromic repeats with the same coordinates, keep one.
    """
    seen = set()
    deduped = []

    for r in repeats:
        if r['type'] == 'F':
            # For forward: normalize so smaller start comes first
            pos1 = (r['qstart'], r['qend'])
            pos2 = (r['sstart'], r['send'])
            key = ('F', tuple(sorted([pos1, pos2])), r['length'])
        elif r['type'] == 'P':
            # For palindromic: normalize coordinates
            pos1 = (min(r['qstart'], r['qend']), max(r['qstart'], r['qend']))
            pos2 = (min(r['sstart'], r['send']), max(r['sstart'], r['send']))
            key = ('P', tuple(sorted([pos1, pos2])), r['length'])
        elif r['type'] == 'R':
            pos1 = (r['qstart'], r['qend'])
            pos2 = (r['sstart'], r['send'])
            key = ('R', tuple(sorted([pos1, pos2])), r['length'])
        elif r['type'] == 'C':
            pos1 = (r['qstart'], r['qend'])
            pos2 = (r['sstart'], r['send'])
            key = ('C', tuple(sorted([pos1, pos2])), r['length'])
        else:
            continue

        if key not in seen:
            seen.add(key)
            deduped.append(r)

    return deduped


def classify_ir_repeats(repeats):
    """Flag large palindromic repeats that are likely IR-derived."""
    for r in repeats:
        if r['type'] == 'P' and r['length'] >= IR_THRESHOLD:
            r['ir_derived'] = True
        else:
            r['ir_derived'] = False
    return repeats


def size_bin(length):
    """Return the size bin label for a given repeat length."""
    for (lo, hi), label in zip(SIZE_BINS, SIZE_BIN_LABELS):
        if lo <= length <= hi:
            return label
    return SIZE_BIN_LABELS[-1]


def process_variety(variety, tmpdir):
    """Process one variety and return results."""
    fasta_path = os.path.join(INPUT_DIR, f"{variety}_std.fasta")
    if not os.path.exists(fasta_path):
        print(f"  WARNING: {fasta_path} not found, skipping.")
        return None

    print(f"  Running BLAST self-alignment for {variety}...")
    repeats, seq_len = run_blast_self(fasta_path, variety, tmpdir)
    print(f"    Raw hits: {len(repeats)}")

    repeats = remove_self_hits(repeats)
    print(f"    After removing self-hits: {len(repeats)}")

    repeats = deduplicate_repeats(repeats)
    print(f"    After deduplication: {len(repeats)}")

    repeats = classify_ir_repeats(repeats)

    return {
        'variety': variety,
        'seq_len': seq_len,
        'repeats': repeats
    }


def write_variety_results(result, output_dir):
    """Write per-variety detailed results to TSV."""
    variety = result['variety']
    repeats = result['repeats']
    outfile = os.path.join(output_dir, f"{variety}_dispersed_repeats.tsv")

    # Sort by type, then by length descending
    repeats_sorted = sorted(repeats, key=lambda r: (r['type'], -r['length']))

    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'Type', 'Length', 'Query_Start', 'Query_End',
            'Subject_Start', 'Subject_End', 'Mismatches', 'Identity%', 'IR_derived'
        ])
        for r in repeats_sorted:
            writer.writerow([
                r['type'], r['length'], r['qstart'], r['qend'],
                r['sstart'], r['send'], r['mismatch'],
                f"{r['pident']:.1f}", 'Yes' if r.get('ir_derived') else 'No'
            ])

    print(f"    Written: {outfile}")


def generate_summary(all_results, output_dir):
    """Generate summary table across all varieties."""
    summary_file = os.path.join(output_dir, "dispersed_repeat_summary.tsv")

    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')

        # ── Table 1: Repeat counts by type ──
        writer.writerow(["=== Dispersed Repeat Counts by Type ==="])
        writer.writerow(["Variety", "Genome_Size", "Total", "Forward(F)", "Palindromic(P)",
                         "Reverse(R)", "Complement(C)", "IR_derived_P"])
        for result in all_results:
            v = result['variety']
            repeats = result['repeats']
            total = len(repeats)
            f_count = sum(1 for r in repeats if r['type'] == 'F')
            p_count = sum(1 for r in repeats if r['type'] == 'P')
            r_count = sum(1 for r in repeats if r['type'] == 'R')
            c_count = sum(1 for r in repeats if r['type'] == 'C')
            ir_count = sum(1 for r in repeats if r.get('ir_derived'))
            writer.writerow([v, result['seq_len'], total, f_count, p_count, r_count, c_count, ir_count])

        writer.writerow([])

        # ── Table 2: Size distribution ──
        writer.writerow(["=== Size Distribution of Dispersed Repeats ==="])
        header = ["Variety"] + SIZE_BIN_LABELS
        writer.writerow(header)
        for result in all_results:
            v = result['variety']
            repeats = result['repeats']
            bins = {label: 0 for label in SIZE_BIN_LABELS}
            for r in repeats:
                bins[size_bin(r['length'])] += 1
            writer.writerow([v] + [bins[label] for label in SIZE_BIN_LABELS])

        writer.writerow([])

        # ── Table 3: Size distribution by type ──
        writer.writerow(["=== Size Distribution by Type (all varieties combined) ==="])
        writer.writerow(["Type"] + SIZE_BIN_LABELS + ["Total"])
        for rtype in ['F', 'P', 'R', 'C']:
            bins = {label: 0 for label in SIZE_BIN_LABELS}
            for result in all_results:
                for r in result['repeats']:
                    if r['type'] == rtype:
                        bins[size_bin(r['length'])] += 1
            total = sum(bins.values())
            type_names = {'F': 'Forward', 'P': 'Palindromic', 'R': 'Reverse', 'C': 'Complement'}
            writer.writerow([f"{type_names[rtype]}({rtype})"] + [bins[label] for label in SIZE_BIN_LABELS] + [total])

        writer.writerow([])

        # ── Table 4: Largest repeats ──
        writer.writerow(["=== Top 10 Largest Repeats (excluding IR-derived) ==="])
        writer.writerow(["Variety", "Type", "Length", "Query_Start", "Query_End",
                         "Subject_Start", "Subject_End", "Mismatches"])
        all_non_ir = []
        for result in all_results:
            for r in result['repeats']:
                if not r.get('ir_derived'):
                    all_non_ir.append((result['variety'], r))
        all_non_ir.sort(key=lambda x: -x[1]['length'])
        for v, r in all_non_ir[:10]:
            writer.writerow([v, r['type'], r['length'], r['qstart'], r['qend'],
                             r['sstart'], r['send'], r['mismatch']])

        writer.writerow([])

        # ── Table 5: IR-derived large repeats ──
        writer.writerow(["=== IR-derived Large Palindromic Repeats ==="])
        writer.writerow(["Variety", "Length", "Query_Start", "Query_End",
                         "Subject_Start", "Subject_End", "Mismatches", "Identity%"])
        for result in all_results:
            for r in result['repeats']:
                if r.get('ir_derived'):
                    writer.writerow([result['variety'], r['length'], r['qstart'], r['qend'],
                                     r['sstart'], r['send'], r['mismatch'], f"{r['pident']:.1f}"])

    print(f"\nSummary written: {summary_file}")
    return summary_file


def print_summary_to_stdout(all_results):
    """Print a concise summary to stdout."""
    print("\n" + "=" * 80)
    print("DISPERSED REPEAT ANALYSIS SUMMARY")
    print(f"Parameters: min_length={MIN_REPEAT_LEN}bp, max_mismatches={MAX_MISMATCHES}")
    print("Types: F=Forward, P=Palindromic, R=Reverse, C=Complement")
    print("=" * 80)

    print(f"\n{'Variety':<8} {'Size':>7} {'Total':>6} {'F':>5} {'P':>5} {'R':>5} {'C':>5} {'IR_P':>5}")
    print("-" * 52)
    for result in all_results:
        repeats = result['repeats']
        f_c = sum(1 for r in repeats if r['type'] == 'F')
        p_c = sum(1 for r in repeats if r['type'] == 'P')
        r_c = sum(1 for r in repeats if r['type'] == 'R')
        c_c = sum(1 for r in repeats if r['type'] == 'C')
        ir_c = sum(1 for r in repeats if r.get('ir_derived'))
        print(f"{result['variety']:<8} {result['seq_len']:>7} {len(repeats):>6} {f_c:>5} {p_c:>5} {r_c:>5} {c_c:>5} {ir_c:>5}")

    print(f"\n{'Size Distribution':}")
    print(f"{'Variety':<8} {'30-44':>7} {'45-59':>7} {'60-89':>7} {'>=90':>7}")
    print("-" * 38)
    for result in all_results:
        bins = defaultdict(int)
        for r in result['repeats']:
            bins[size_bin(r['length'])] += 1
        print(f"{result['variety']:<8} {bins['30-44']:>7} {bins['45-59']:>7} {bins['60-89']:>7} {bins['>=90']:>7}")

    # Largest non-IR repeats
    print(f"\nTop 5 largest repeats (excluding IR-derived):")
    all_non_ir = []
    for result in all_results:
        for r in result['repeats']:
            if not r.get('ir_derived'):
                all_non_ir.append((result['variety'], r))
    all_non_ir.sort(key=lambda x: -x[1]['length'])
    for v, r in all_non_ir[:5]:
        print(f"  {v}: {r['type']} {r['length']}bp at {r['qstart']}-{r['qend']} / {r['sstart']}-{r['send']} ({r['mismatch']} mismatches)")

    # IR repeats
    print(f"\nIR-derived large palindromic repeats:")
    for result in all_results:
        for r in result['repeats']:
            if r.get('ir_derived'):
                print(f"  {result['variety']}: P {r['length']}bp at {r['qstart']}-{r['qend']} / {r['sstart']}-{r['send']}")


def main():
    print("Dispersed Repeat Analysis (REPuter-equivalent via BLASTn)")
    print(f"Parameters: min_length={MIN_REPEAT_LEN}bp, max_mismatches={MAX_MISMATCHES}")
    print(f"Input: {INPUT_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for variety in VARIETIES:
            print(f"Processing {variety}...")
            result = process_variety(variety, tmpdir)
            if result:
                write_variety_results(result, OUTPUT_DIR)
                all_results.append(result)
                print(f"    Total dispersed repeats: {len(result['repeats'])}")
            print()

    if all_results:
        generate_summary(all_results, OUTPUT_DIR)
        print_summary_to_stdout(all_results)

    print("\nDone!")


if __name__ == "__main__":
    main()
