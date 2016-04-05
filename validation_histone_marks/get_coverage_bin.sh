#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt
for dir in /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_input/*/
do
    dir=${dir%*/}
    #echo $dir
    dir_only=${dir##*/}
    mkdir "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only

    for dir2 in $dir/*
    do
  
	dir2_only=${dir2##*/}
	mkdir "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only"/"$dir2_only

	peakFile="/data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir2_only"/peaks/filtered.peaks.narrowPeak"
	echo $peakFile
	
	
	#loop on all mid files 
	for DHS_file in /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/*
	do
		DHS_file_only=${DHS_file##*/}
		
		outputFile="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only"/"$dir2_only"/"$DHS_file_only	
		#echo $DHS_file
		#echo $outputFile		
		perl BedToPeakCoverageMeta.pl $peakFile $DHS_file $outputFile 300
	done
#done
	
done
done


