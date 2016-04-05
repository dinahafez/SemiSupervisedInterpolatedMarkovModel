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
if [ $dir_only != "H3K27ac" ] 
then
    Dir_output="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only
    cmd="cat "

    for dir2 in $dir/*
    do
  
	dir2_only=${dir2##*/}
	peakFile="/data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir2_only"/peaks/all.peaks.narrowPeak"
	cmd=$cmd$peakFile" ";
    done

    cmd=$cmd" > /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir_only".all.peaks.NarrowPeaks"		
#echo $cmd


#loop on all mid files 
for DHS_file in /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/*
do
	DHS_file_only=${DHS_file##*/}
		
	outputFile="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only"/"$DHS_file_only".filtered"
	peakFile="/data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir_only".filtered.peaks.NarrowPeaks"
	echo $outputFile
		#echo $DHS_file
		#echo $outputFile		
		perl BedToPeakCoverageMeta.pl $peakFile $DHS_file $outputFile 10
	done



fi




	
done


