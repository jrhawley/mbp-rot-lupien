#!/bin/bash

#$ -N MuSE
#$ -cwd 
#$ -q lupiengroup
module load perl/5.18.1  
module load R

#mutations.sorted.bed C3D.output DHS_reference.bed window cor
perl MuSE.pl -m $1 -c $2 -b $3 -w $4 -t $5 
