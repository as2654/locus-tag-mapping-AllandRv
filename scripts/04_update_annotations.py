#!/usr/bin/env python3
import pandas as pd
import re
from tqdm import tqdm
import argparse
import sys
import os

# ---------------------------
# Step 1: Parse arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Map new locus tags to old ones using BLAST and update GTF/GFF files.")
parser.add_argument("--blast_tab", required=True, help="BLAST output .tab file (outfmt 6 with qcovs, scovs, slen)")
parser.add_argument("--gtf_file", required=False, help="GTF file to update locus tags")
parser.add_argument("--gff_file", required=False, help="GFF file to update locus tags")
parser.add_argument("--output_dir", default=".", help="Directory to write output mapping and updated files (default is current directory)")
parser.add_argument("--min_pident", type=float, default=70.0, help="Minimum percent identity to consider a hit (default=70)")
parser.add_argument("--min_qcovs", type=float, default=90.0, help="Minimum query coverage to consider a hit (default=90)")
parser.add_argument("--min_scovs", type=float, default=90.0, help="Minimum subject coverage to consider a hit (default=90)")
args = parser.parse_args()

blast_tab = args.blast_tab
gtf_file = args.gtf_file
gff_file = args.gff_file
output_dir = args.output_dir
min_pident = args.min_pident
min_qcovs = args.min_qcovs
min_scovs = args.min_scovs

os.makedirs(output_dir, exist_ok=True)

# ---------------------------
# Step 2: Read BLAST results
# ---------------------------
print("Reading BLAST results...")
blast_cols = ["query", "subject", "pident", "qcovs", "scovs", "length", "slen", "evalue", "bitscore"]
blast_results = pd.read_csv(blast_tab, sep='\t', names=blast_cols)

# ---------------------------
# Step 3: Filter BLAST hits by thresholds
# ---------------------------
print(f"Filtering BLAST hits: pident>={min_pident}, qcovs>={min_qcovs}, scovs>={min_scovs}")
filtered_results = blast_results[
    (blast_results.pident >= min_pident) &
    (blast_results.qcovs >= min_qcovs) &
    (blast_results.scovs >= min_scovs)
]

# ---------------------------
# Step 4: Determine best hit per query
# ---------------------------
print("Selecting best hits per query...")
new_to_old_mapping = {}
for query, group in filtered_results.groupby('query'):
    best_hit = group.sort_values('bitscore', ascending=False).iloc[0]
    new_to_old_mapping[query] = best_hit['subject']

# ---------------------------
# Step 5: Number duplicates
# ---------------------------
print("Handling duplicates in old locus tags...")
numbered_mapping = {}
old_tag_counts = {}
for new_tag, old_tag in new_to_old_mapping.items():
    old_tag_counts[old_tag] = old_tag_counts.get(old_tag, 0) + 1

old_tag_current_number = {}
for new_tag, old_tag in new_to_old_mapping.items():
    count = old_tag_counts[old_tag]
    if count > 1:
        old_tag_current_number[old_tag] = old_tag_current_number.get(old_tag, 0) + 1
        numbered_old_tag = f"{old_tag}.{old_tag_current_number[old_tag]}"
    else:
        numbered_old_tag = old_tag
    numbered_mapping[new_tag] = numbered_old_tag

# ---------------------------
# Step 6: Save mapping TSV
# ---------------------------
mapping_file = os.path.join(output_dir, "locus_tag_mapping.tsv")
print(f"Writing mapping to {mapping_file}")
with open(mapping_file, "w") as f:
    f.write("new_locus_tag\told_locus_tag\n")
    for new_tag, numbered_old_tag in numbered_mapping.items():
        f.write(f"{new_tag}\t{numbered_old_tag}\n")

# ---------------------------
# Step 7: Functions to replace locus_tag
# ---------------------------
def replace_locus_tag(attr_str, mapping):
    pattern = r'locus_tag[= ]"?(.*?)"?(;|$)'
    def repl(match):
        old_value = match.group(1)
        sep = match.group(2)
        new_value = mapping.get(old_value, old_value)
        return f'locus_tag={new_value}{sep}'
    return re.sub(pattern, repl, attr_str)

def replace_locus_tag_gtf(attr_str, mapping):
    pattern = r'locus_tag "([^"]+)"'
    def replacer(match):
        old_value = match.group(1)
        return f'locus_tag "{mapping.get(old_value, old_value)}"'
    return re.sub(pattern, replacer, attr_str)

# ---------------------------
# Step 8: Process GTF file
# ---------------------------
if gtf_file:
    print(f"Updating GTF file: {gtf_file}")
    gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#')
    gtf_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 
                      'score', 'strand', 'frame', 'attribute']
    gtf_df['attribute'] = gtf_df['attribute'].apply(lambda x: replace_locus_tag_gtf(x, numbered_mapping))
    gtf_outfile = os.path.join(output_dir, os.path.basename(gtf_file).replace(".gtf", "_updated.gtf"))
    gtf_df.to_csv(gtf_outfile, sep='\t', header=False, index=False, quoting=3)
    print(f"Updated GTF written to: {gtf_outfile}")

# ---------------------------
# Step 9: Process GFF file
# ---------------------------
if gff_file:
    print(f"Updating GFF file: {gff_file}")
    gff_lines = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                gff_lines.append(line)
            else:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 9:
                    parts[8] = replace_locus_tag(parts[8], numbered_mapping)
                    gff_lines.append("\t".join(parts) + "\n")
                else:
                    gff_lines.append(line)
    gff_outfile = os.path.join(output_dir, os.path.basename(gff_file).replace(".gff", "_updated.gff"))
    with open(gff_outfile, "w") as f:
        f.writelines(gff_lines)
    print(f"Updated GFF written to: {gff_outfile}")

# ---------------------------
# Step 10: Summary
# ---------------------------
print(f"Total new locus tags mapped: {len(new_to_old_mapping)}")
print(f"Unique old locus tags used: {len(set(new_to_old_mapping.values()))}")
print(f"Old locus tags with duplications: {sum(1 for tag,count in old_tag_counts.items() if count > 1)}")


