#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o errorR.out -j y
#$ -e errorR.log

Rscript /data/ohler/Scott/R_Scripts/Meta.R /data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/H3K27ac/11_initialize_redfly_mid.bed.filtered.peakCoverageMeta

