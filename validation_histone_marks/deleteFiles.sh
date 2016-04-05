#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt
for dir in /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/*/
do
    dir=${dir%*/}   #Histone mark name

    echo $dir
    dir_only=${dir##*/}

    Dir_output="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$dir_only	
   rm $dir/*initialize*
  	
	#mkdir $Dir_output"/Overlap_Bin10/"
	#cmd="mv "$Dir_output"/*.tiff "$Dir_output"/Overlap_Bin10/"
	#$cmd
done
