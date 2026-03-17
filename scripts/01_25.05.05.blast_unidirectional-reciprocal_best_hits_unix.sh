#!/bin/bash

# Avi Shah - W. Evan Johnson Lab
# The purpose of this script is to perform a unidirectional best hits 
# analysis. This is edited from an older reciprocal best hits script.
# We will make a blast database of the old H37Rv reference, and find 
# the best hits unidirectionally from the old reference mapping to each 
# gene feature in the new Alland reference to eventually make the 
# appropriate replacements for the locus tag features in the Alland 
# reference's GTF file.

# Check for two input parameters
if [ ${#} -ne 2 ]; then
    echo "Usage: ${0} <new_reference.fasta> <old_reference.fasta>"
    exit 1
fi

new_fasta=${1}
old_fasta=${2}

# Move to blast_resources directory (relative to script location)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}/../blast_resources" || exit 1

# Create database prefixes by removing .fasta extension
new_db="$(basename "${new_fasta%.*}")"
old_db="$(basename "${old_fasta%.*}")"

# Create BLAST database
makeblastdb -in "${old_fasta}" -dbtype nucl -out "${old_db}" || exit 1

# Run unidirectional BLAST searches with stringent criteria
blastn -query "${new_fasta}" -db "${old_db}" \
  -outfmt "6 qseqid sseqid pident qcovs length slen evalue bitscore" \
  -perc_identity 25 \
  -qcov_hsp_perc 50 \
  -evalue 1e-2 \
  -out "${new_db}_vs_${old_db}.tab" \
  -num_threads $(nproc) || exit 1

# Add subject coverage statistic as 5th column
tmp=$(mktemp)
awk 'BEGIN{OFS="\t"} {
    scovs = ($5 / $6) * 100;
    print $1, $2, $3, $4, scovs, $5, $6, $7, $8
}' "${new_db}_vs_${old_db}.tab" > "$tmp" && mv "$tmp" "${new_db}_vs_${old_db}.tab"

echo "BLAST results saved to:"
echo "- ${new_db}_vs_${old_db}.tab"
