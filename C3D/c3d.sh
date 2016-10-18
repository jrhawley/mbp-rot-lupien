#!/bin/bash
DIR=/mnt/work1/users/lupiengroup/TahmidProject/C3D/
#$ -cwd
#$ -N C3D
#$ -o logs/c3d.log
#$ -q lupiengroup
#$ -S /bin/bash

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi
# Princess Margaret Cancer Centre - University Health Network, July 25, 2016

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Takes correlations between open regions of chromatin based on DNaseI hypersensitivity signals
# Regions with high correlations are candidates for 3D interactions
# Performs association tests on each candidate & adjusts p-values
# Identifies transcription factor motifs which overlap open regions
# Produces interaction landscapes and motif tracks in PDF format

# an array for parameters
typeset -A config
# default values
config=(
    	[reference]=""
    	[db]=""
    	[anchor]=""
	[outDirectory]=""
	[matrix]=""
    	[window]="500000"
    	[correlationThreshold]="0.5"
	[pValueThreshold]="0.05"
    	[qValueThreshold]="0.05"
	[correlationMethod]="pearson"
    	[motifs]=""
	[fimoPValue]="0.0001"
    	[fimoQValue]="1"
    	[figures]=""
	[figureWidth]="500000"
	[zoom]="0"
	[colours]="#bdd7e7,#6baed6,#3182bd,#08519c"
	[tracks]="y"
	[sampleName]=" "
)
# read parameters from config file
while read line
do
	if echo $line | grep -F = &>/dev/null
    	then
		eval expanded_line="$line"
        	varname=$(echo "$line" | cut -d '=' -f 1)
        	config[$varname]=$(echo $expanded_line | cut -d '=' -f 2-)
    	fi
	if echo $line | grep -F 'module load' &>/dev/null
    	then
    		eval $line
	fi
done < $1

# if one of these parameters is empty, they'll get their default values
if [ -z "${config[window]}" ]; then
        config[window]="500000"
fi
if [ -z "${config[correlationThreshold]}" ]; then
        config[correlationThreshold]="0.5"
fi
if [ -z "${config[pValueThreshold]}" ]; then
        config[pValueThreshold]="0.05"
fi
if [ -z "${config[qValueThreshold]}" ]; then
        config[qValueThreshold]="0.05"
fi
if [ -z "${config[correlationMethod]}" ]; then
        config[correlationMethod]="pearson"
fi
if [ -z "${config[fimoPValue]}" ]; then
        config[fimoPValue]="0.0001"
fi
if [ -z "${config[fimoQValue]}" ]; then
        config[fimoQValue]="1"
fi
if [ -z "${config[figureWidth]}" ]; then
        config[figureWidth]=${config[window]}
fi
if [ -z "${config[zoom]}" ]; then
        config[zoom]="0"
fi
if [ -z "${config[colours]}" ]; then
        config[colours]="#bdd7e7,#6baed6,#3182bd,#08519c"
fi
if [ -z "${config[tracks]}" ]; then
        config[tracks]="y"
fi
if [ -z "${config[sampleName]}" ]; then
        config[sampleName]=" "
fi
# the index of the sample
trackNumber=1
# the number of total samples
numSamples=1
# overwrite parameters for multiple samples
if [ "$2" = "-ref" ]; then
        config[reference]="$3"
        config[matrix]=""
fi
if [ "$2" = "-matrix" ]; then
        config[matrix]="$3"
fi
if [ "$4" = "-out" ]; then
        config[outDirectory]="$5"
fi
if [ "$6" = "-sample" ]; then
        config[sampleName]="$7"
fi
if [ "$8" = "-track" ]; then
        trackNumber="$9"
fi
if [ "${10}" = "-numSamples" ]; then
        numSamples="${11}"
fi

# make output directory
mkdir -p ${config[outDirectory]}

# timestamp function
timestamp() {
date +"%Y-%m-%d_%H-%M-%S"
}

# if the anchor file does not have 5 fields, create a new formatted one
anchorCols=$(awk '{print NF}' ${config[anchor]} | sort -nu | tail -n 1)
if [ "$anchorCols" -ne "5" ]; then
        awk '{print $1"\t"$2"\t"$3"\t.\t"$1":"$2"-"$3}' ${config[anchor]} > ${config[outDirectory]}/anchors.bed
else
	cat ${config[anchor]} > ${config[outDirectory]}/anchors.bed
fi
config[anchor]="${config[outDirectory]}/anchors.bed"

# if matrix is missing, map bedgraphs to reference
if [ "${config[matrix]}" = "" ]; then
        # Map all background files to one reference sample/peak catalogue
        echo "$(timestamp): Mapping peak files"
        cut -f 1-3 ${config[reference]} > ${config[outDirectory]}/ref.bed

        counter=1
        cat ${config[db]} | while read i; do
                echo "Mapping ${i} to file $counter.map.bed"
                mapBed -a ${config[outDirectory]}/ref.bed -b ${i} -c 4 -o max -null 0 | awk 'BEGIN{ OFS="\t" }{ print $1, $2, $3, $4 }' > ${config[outDirectory]}/$counter.map.bed
                counter=$((counter + 1))
        done
else
	tail -n +2 ${config[matrix]} | awk '{print $1}' | awk -F'[:-]' '{print $1"\t"$2"\t"$3}' > ${config[outDirectory]}/ref1.bed
	sort -k1,1 -k2,2n ${config[outDirectory]}/ref1.bed > ${config[outDirectory]}/ref.bed
	rm ${config[outDirectory]}/ref1.bed
fi

# Run R script
Rscript $DIR/c3d.R "${config[outDirectory]}" "${config[outDirectory]}" "${config[anchor]}" "${config[db]}" "${config[window]}" "${config[correlationThreshold]}" "${config[pValueThreshold]}" "${config[qValueThreshold]}" "${config[correlationMethod]}" "${config[matrix]}" "${config[motifs]}" "${config[fimoPValue]}" "${config[fimoQValue]}" "${config[figures]}" "${config[figureWidth]}" "${config[zoom]}" "${config[colours]}" "${config[tracks]}" "${config[sampleName]}" "$trackNumber" "$numSamples" "$(timestamp)" "$DIR"

