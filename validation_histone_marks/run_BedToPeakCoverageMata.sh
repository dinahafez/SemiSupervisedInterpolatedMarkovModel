#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error.out -j y
#$ -l h_vmem=20G
#$ -l mem_free=20G
#$ -e error.log

perl BedToPeakCoverageMeta.pl $1 $2 $3 $4
