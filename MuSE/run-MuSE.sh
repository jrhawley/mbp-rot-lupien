#!/bin/bash

#$ -N MuSE
#$ -cwd 
#$ -q lupiengroup
module load perl/5.18.1  
module load R

#mutations.sorted.bed C3D.output DHS_reference.bed window cor
perl MuSE-1.6.pl $1 $2 $3 $4 $5 
