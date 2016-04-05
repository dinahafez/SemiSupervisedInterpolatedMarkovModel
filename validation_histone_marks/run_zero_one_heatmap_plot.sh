for dir in /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/*/
do
    dir=${dir%*/}

    #echo $dir
    dir_only=${dir##*/}
	echo $dir 

    for coverage_file in $dir/*.filtered.peakCoverageMeta
    do
	#	coverage_file_only=${coverage_file##*/}
	#echo $coverage_file
#	Rscript /data/ohler/Dina/Drosophila/Research/code/validaiton_histone_marks/Zero_One_HeatMap.R $coverage_file\
	Rscript /data/ohler/Dina/Drosophila/Research/code/validaiton_histone_marks/MetaLine.R $coverage_file
    done
done


