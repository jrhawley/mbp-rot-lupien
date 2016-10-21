#!/bin/bash
#$ -cwd
#$ -N makeanchors
#$ -j y
#$ -q lupiengroup
#$ -S /bin/bash


#1 = list of genes, 1 gene per line
#2 = number of bases upstream of TSS
#3 = number of bases downstream of TSS
#4 = output file

module load R/3.2.2

Rscript make_anchors_from_genes.R "$1" "$2" "$3" "$4"


