#!/bin/bash

#$ -S /bin/bash
#$ -N MuSE
#$ -j y
#$ -cwd 
#$ -q light.q 
module load perl/5.18.1  
module load R

#mutations.sorted.bed IFC.output DHS.bed window cor

perl MUSE-BETA-1.5-1-mutOnly.pl $1 $2 $3 $4 $5 
