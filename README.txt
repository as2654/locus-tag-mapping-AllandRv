We need to replace locus tags in the new reference of the Mycobacterium tuberculosis
H37Rv reference genome because existing functional/pathway analysis tools rely on 
'Rv numbers' which refer to originally sequenced H37Rv strain from 1998, and that was 
updated in 2013.

I will use the 2013 update to the H37Rv reference to modify the locus tags of the 
updated H37Rv reference genome here.

You will need GenBank v3 (.gbff) formatted annotation/genome files for each reference
to run this pipeline.

You should run the executable notebooks/scripts in the order of the prefixes (0X_).


The first script extracts each gene features nucleotide sequence (with appropriate 
reverse complementation) and arranges them into a multi-FASTA file.


The second script makes a blast database using the old reference and runs a blastn 
search for each gene in the new reference against the old, finding the best hits 
for a given percent identity and query cover threshold (originally 97% for both, 
25% percent identity and 50% query cover in the current implementation as we want 
to map as many locus tags as comprehensively as possible).

Briefly, 3817 of the 4164 locus tags extracted from the Alland reference had matching 
nucleotide sequences in the old reference, with mappings to 3792 unique old locus 
tags, when performing a stringent 97/97 search. There were 23 unique old sequences 
that multi-mapped in this manner. All multi-mapping annotations are saved in a 
stringent_multi-mappers.tsv file.

After loosening the requirements to 25% identity and 50% query cover with a weaker 
e-value threshold (1e-2 versus the default 1e-5), we found mappings for 3991 new 
genes from 3945 old genes, with 39 multi-mapping old locus tags. These same hits 
were retained when increasing stringency to 30% identity and 75% query cover.

An intermediate stringency was also attempted (70% identity and 90% query cover at 
an e-value threshold of 1e-2) which yielded 3898 replacements from 3864 old locus 
tags and 31 multi-mappers.


The third script performs a best hits analysis, generating a tab-separated value 
file with two columns. The first column has each gene feature's locus tag from the 
Alland reference, and the second column has the best replacement locus tag from 
the old reference, defaulting to the same locus tag from the Alland reference if a 
hit with our minimum percent identity and query cover requirements is not found.
The actual process of identifying best hits has been simplified to simply using 
the hit with the highest bitscore (last column of blast .tab results file), as it 
is a measure of statistical significance of the alignment, and better than trying 
to use one of the values (or some ad-hoc composite of the various scores).
When the same gene in the old reference is the best match to multiple new genes, 
its locus tag will have a counter integer appended as a suffix (.N) so that each 
locus tag will be unique in both columns.


The fourth script parses this TSV file and replaces the locus tag attribute for 
each feature in the last column of the GTF annotations file for the new reference.
It will do so comprehensively (for gene, transcript, etc. features), but at least 
for each gene feature, which should make the outputs for counting workflows (using 
e.g. featureCounts) compatible with old pathway tools.
