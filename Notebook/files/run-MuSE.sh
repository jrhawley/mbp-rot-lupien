#!/bin/bash

#$ -N MuSE
#$ -cwd 
#$ -q lupiengroup
#$ -t 1-8
module load perl/5.18.1  
module load R

GENES=(GATA3 FOXA1 ESR1 PGR ERBB2 BRCA1 BRCA2 TP53)
IDX=$((${SGE_TASK_ID}-1))
C3D_FILE="$2.${GENES[$IDX]}"


#mutations.sorted.bed C3D.output DHS_reference.bed window cor
perl MuSE.pl -m $1 -c ${C3D_FILE} -b $3 -w $4 -t $5 
