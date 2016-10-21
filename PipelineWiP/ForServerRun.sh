#!/bin/bash
#$ -cwd
#$ -V
#$ -j yes
#$ -pe smp 8
python forServerRun.py -G Phyca11_filtered_genes.gff > runInfo.txt