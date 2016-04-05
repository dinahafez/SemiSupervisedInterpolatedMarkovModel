#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt

unlabeledDir="/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/"
posDir="/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output/"
outputDir="/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output/TF_enrichment/";

for i in {11..39..1}
do

	intersectBed -a $unlabeledDir$i"_unlabeled_in_uniq.bed" -b $posDir$i"_pos_stage_2.bed" -v > $outputDir$i"unlabled_not_pos.bed";
	intersectBed -a $unlabeledDir$i"_specified_redfly_uniq.bed" -b $posDir$i"_pos_stage_2.bed" -v > $outputDir$i"_initialize_not_pos.bed";
	cat $outputDir$i"unlabled_not_pos.bed"  $outputDir$i"_initialize_not_pos.bed" >  $outputDir$i"_in_50k_not_pos.bed";
done
