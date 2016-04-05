#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt
for dir in /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_input/*/
do
    dir=${dir%*/}   #Histone mark name

    #echo $dir
    dir_only=${dir##*/}

    #Dir_output="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only	
    cmd="cat "

    for dir2 in $dir/*
    do
  	
	dir2_only=${dir2##*/}
	if [ $dir2_only != "Embryo_16_20hr" ] &&  [ $dir2_only != "Embryo_0_4hr" ];
	then	
	peakFile="/data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir2_only"/peaks/filtered.peaks.narrowPeak"
	cmd=$cmd$peakFile" ";
	fi
   done

    cmd=$cmd" > /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir_only".4.12.filtered.peaks.NarrowPeaks"
#echo $cmd
#eval "$cmd" 




#loop on all mid files 
#for DHS_file in /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/*.bed
#do
	DHS_file="/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_no_overlap_TSS_not_pos_14_mid.bed"

	#"/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/DHS_no_overlap_TSS_not_common_12_mid.bed"

	DHS_file_only=${DHS_file##*/}
	mkdir "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/"$dir_only

	outputFile="/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/"$dir_only"/"$DHS_file_only".filtered"
	peakFile="/data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/"$dir_only"/"$dir_only".4.12.filtered.peaks.NarrowPeaks"  #all peaks for one histone mark across all timepoits
	#strandFile="/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed"
	#echo $outputFile
	echo $DHS_file
	#echo $outputFile		
	#if [[ $DHS_file_only == *pos* ]]
	#then
		qsub run_BedToPeakCoverageMata.sh $peakFile $DHS_file $outputFile 10
#		qsub run_BedToPeakCoverageMata.sh $peakFile $DHS_file $outputFile 10 $strandFile
	#fi	
	done	
#done
