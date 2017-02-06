#!/bin/bash
set -e



# We got the original bins in 00_BINS_FASTA from here: http://ggkbase.berkeley.edu/carrol,
# and stored them in individual FASTA files and archived the entire directory as `01_BINS_FASTA.tar.gz`.
# following steps identifies the annotation for each split in Sharon et al results, and creates an
# 'additional metadata file' that can be imported from the interactive interface to visualize the
# concordance.
export contigs_db="/PATH/TO/CONTIGS.db"

# Clean stuff from the previous analysis:
    ./00_CLEAN.sh

# Open the archive

    tar -zxvf 01_BINS_FASTA.tar.gz

# Create one FASTA file for all Sharon et al bins:

    cd 01_BINS_FASTA/ && ls *.fa | awk 'BEGIN{FS="."}{print $1}' > ../02_SHARON_et_al_BIN_IDS.txt && cd ..
    rm -rf 03_SHARON_et_al_BINS_FASTA.fa
    for b in `cat 02_SHARON_et_al_BIN_IDS.txt`
    do
        echo '>'$b >> 03_SHARON_et_al_BINS_FASTA.fa
        grep -v '>' 01_BINS_FASTA/$b.fa >> 03_SHARON_et_al_BINS_FASTA.fa
    done

# Run this to get contigs

    anvi-export-contigs -c $contigs_db -o 04_INFANT_CONTIGS.fa

# Create a BLAST db for Sharon et al bins:

    makeblastdb -in 03_SHARON_et_al_BINS_FASTA.fa -dbtype nucl

# Do the BLAST search:

    blastn -query 04_INFANT_CONTIGS.fa -db 03_SHARON_et_al_BINS_FASTA.fa -out 05_SHARON_et_al_BINS_HITS.b6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -max_target_seqs 1

# Create an additional metadata file for anvi'o:

    echo "contigs	Sharon_et_al_bin" > 06_ADDITIONAL_METADATA_FOR_SHARON_et_al_BINS.txt
    awk '{if($3 > 95) print $1"	"$2}' 05_SHARON_et_al_BINS_HITS.b6 | sort | uniq | sort -k 2 >> 06_ADDITIONAL_METADATA_FOR_SHARON_et_al_BINS.txt
