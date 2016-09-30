#!/bin/bash
#$ -cwd
#$ -V
#$ -j yes
#$ -pe smp 8
python forServerRun.py -Q Y006W_S7_L001_R1_001P100.fastq -Q Y006W_S7_L001_R2_001P100.fastq -Q AD84_S8_L001_R1_001P100.fastq -Q AD84_S8_L001_R2_001P100.fastq -Q Q108_S9_L001_R1_001P100.fastq -Q Q108_S9_L001_R2_001P100.fastq -g Phyca11_unmasked_genomic_scaffolds.fasta -g LGSJ01.1.fsa_nt -g LGTR01.1.fsa_nt -g AUPN01.1.fsa_nt -g AUVH01.1.fsa_nt -G Phyca11_filtered_genes.gff > runInfo.txt