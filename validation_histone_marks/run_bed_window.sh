for dir in /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/*/
do
    dir=${dir%*/}
    dir_only=${dir##*/}
    mkdir /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/$dir_only/overlap_count/

    for chip_file in $dir/*.all.peaks.NarrowPeaks
    do
	#	coverage_file_only=${coverage_file##*/}
	#echo $coverage_file

	DHS_file="/data/ohler/Dina/Drosophila/Research/data/DHS/DHS_overlap_TSS_10bp.bed"
#DHS_file="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/ubiquitous.bed"
#	for DHS_file in /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/MC_5order/MC_ubiq_output/*.bed
#	do
	      DHS_file_only=${DHS_file##*/}
              outputFile="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only"/overlap_count/"$DHS_file_only
               bedtools window -w 75 -u -a $DHS_file -b $chip_file > $outputFile

#	        done
    done
done

