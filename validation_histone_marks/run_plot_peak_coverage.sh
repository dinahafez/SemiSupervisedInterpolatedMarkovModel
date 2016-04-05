for dir in /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/*/
do
    dir=${dir%*/}

    #echo $dir
    dir_only=${dir##*/}
	echo $dir 

    for coverage_file in $dir/*.peakCoverageMeta
    do
	#	coverage_file_only=${coverage_file##*/}
	#echo $coverage_file
	Rscript /data/ohler/Dina/Drosophila/Research/code/validaiton_histone_marks/Meta.R $coverage_file
    done
done


