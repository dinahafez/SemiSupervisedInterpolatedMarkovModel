#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt

#convert Fasta to bed
#bedmap
#sum of average 

inputDir=$1;
outputDir="/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output/TF_enrichment/";

for TF in /data/ohler/Dina/Drosophila/Research/data/TF_chip/*.bed
do 
	TF=${TF%*/}
	TF_only=${TF##*/}
	#echo $TF_only

	#for file in $inputDir*_pos_stage_2.bed
	for i in {11..39..1}
	do
		file=$inputDir$i"_pos_stage_2.bed"   #$outputDir$i"_initialize_not_pos.bed"   #"unlabled_not_pos.bed"
		file_only=${file##*/}
		#echo $file
		if [ ! -s $outputDir$TF_only"_"$file_only ]
		then
			bedmap --delim "\t" --faster --mean $file $TF > $outputDir$TF_only"_"$file_only
		fi
		awk '{if ( $1 >=0 ) s+=$1   } END {print s}' $outputDir$TF_only"_"$file_only
		#wc -l  $file
	done
done
