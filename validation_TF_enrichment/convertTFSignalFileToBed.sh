#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt
for dir in /data/ohler/Dina/Drosophila/Research/data/TF_chip/*.wig
do
    dir=${dir%*/}
	   

 	echo $dir
   	dir_only=${dir##*/}
	wig2bed < $dir > $dir".bed"
	
done
