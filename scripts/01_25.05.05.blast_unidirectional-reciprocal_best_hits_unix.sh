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

# Create database prefixes by removing .fasta extension
new_db="${new_fasta%.*}"
old_db="${old_fasta%.*}"

# Create BLAST databases
#makeblastdb -in "${new_fasta}" -dbtype nucl -out "${new_db}" || exit 1
makeblastdb -in "${old_fasta}" -dbtype nucl -out "${old_db}" || exit 1

# Run unidirectional BLAST searches with stringent criteria
blastn -query "${new_fasta}" -db "${old_db}" \
  -outfmt "6 qseqid sseqid pident qcovs length evalue bitscore" \
  -perc_identity 70 \
  -qcov_hsp_perc 90 \
  -evalue 1e-2 \
  -out "${new_db}_vs_${old_db}.tab" \
  -num_threads $(nproc) || exit 1

#blastn -query "${old_fasta}" -db "${new_db}" \
#  -outfmt "6 qseqid sseqid pident qcovs length evalue bitscore" \
#  -perc_identity 97 \
#  -qcov_hsp_perc 97 \
#  -out "${old_db}_vs_${new_db}.tab" \
#  -num_threads $(nproc) || exit 1

echo "BLAST results saved to:"
echo "- ${new_db}_vs_${old_db}.tab"
#echo "- ${old_db}_vs_${new_db}.tab"
