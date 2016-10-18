#!/bin/bash
#$ -cwd
#$ -N makeTracks
#$ -o logs/makeTracks.log
#$ -q lupiengroup
#$ -S /bin/bash

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi
# Princess Margaret Cancer Centre - University Health Network, July 22, 2016

# Creates a custom track by combining outputs of C3D for a single sample or multiple samples & appends motif occurances. This track can be viewed on the genome browser.
module load bedtools
# parameters
anchor="$1"
outDirectory="$2"
references="$3"
motifs="$4"
# an array of genes from the last column of the anchor file
genes=( $(cut -d$'\t' -f5 $anchor ) )
# if references is not passed, then there's only 1 sample
if [ "$references" = "" ]; then
	# iterate through genes & make custom track for each gene
        for i in "${!genes[@]}"; do
        	cat $outDirectory/${genes[$i]}.anchor $outDirectory/${genes[$i]}.bedGraph $outDirectory/M*.bed > $outDirectory/${genes[$i]}.tracks.txt
	done
else # multiple samples
	# an array of bed files for each sample
        files=( $(cut -d ' ' -f1 $references ) )
	# an array of sample names
	sample=( $(cut -d ' ' -f2 $references ) )
	# if motifs were provided
	if [ "$motifs" != "" ]; then
		# an array of JASPAR IDs for the TFs in the motifs file
		motifIDs=( $(cut -d ' ' -f1 $motifs ) )
		for i in "${!motifIDs[@]}"; do
        		for j in "${!sample[@]}"; do
                		motifFiles[$j]="$outDirectory/${sample[$j]}/${motifIDs[$i]}.bed"
        		done
			# concatenate all motif BEDs from each sample
        		cat ${motifFiles[*]} > $outDirectory/${motifIDs[$i]}.1.txt
			# remove duplicate coordinates
			awk '!a[$0]++' $outDirectory/${motifIDs[$i]}.1.txt > $outDirectory/${motifIDs[$i]}.2.txt
			# sort the BED
        		(head -n 1 $outDirectory/${motifIDs[$i]}.2.txt && tail -n +2 $outDirectory/${motifIDs[$i]}.2.txt | sort -k1,1 -k2,2n) > $outDirectory/${motifIDs[$i]}.bed
			rm $outDirectory/${motifIDs[$i]}.1.txt $outDirectory/${motifIDs[$i]}.2.txt
		done
	fi
	
	for i in "${!sample[@]}"; do
        	refFiles[$i]="$outDirectory/${sample[$i]}/ref.bed"
        done
	bedtools multiinter -i ${refFiles[*]} | awk -v numSamples="${#files[@]}" '$4 == numSamples' | awk '{print $1"\t"$2"\t"$3}' > $outDirectory/union.bed
	echo 'track name="Sample Intersections" description=" " visibility=1 color=0,0,0 db=hg19' | cat - $outDirectory/union.bed > $outDirectory/temp.txt && mv $outDirectory/temp.txt $outDirectory/union.bed
	
	for i in "${!genes[@]}"; do
		# for each gene, write the browser header & anchor coordinates to a file
        	cat $outDirectory/${sample[0]}/${genes[$i]}.anchor > $outDirectory/${genes[$i]}.1.txt
        	for j in "${!sample[@]}"; do
                	geneFiles[$j]="$outDirectory/${sample[$j]}/${genes[$i]}.bedGraph"
        	done
		# concatenate the BEDGRAPHs for each sample & the union.bed to the anchor file
        	cat $outDirectory/${genes[$i]}.1.txt ${geneFiles[*]} $outDirectory/union.bed > $outDirectory/${genes[$i]}.2.txt
		
		if [ "$motifs" != "" ]; then
        		for j in "${!motifIDs[@]}"; do
                		motifFiles[$j]="$outDirectory/${motifIDs[$j]}.bed"
        		done
			# append motif coordinates of each TF to the file
        		cat $outDirectory/${genes[$i]}.2.txt ${motifFiles[*]} > $outDirectory/${genes[$i]}.tracks.txt
			rm $outDirectory/${genes[$i]}.2.txt $outDirectory/${genes[$i]}.1.txt
		else # if motifs weren't provided
			# just copy the concatenated BEDGRAPHs to a better named file
			cp $outDirectory/${genes[$i]}.2.txt $outDirectory/${genes[$i]}.tracks.txt
			rm $outDirectory/${genes[$i]}.2.txt $outDirectory/${genes[$i]}.1.txt
		fi
	done
fi

