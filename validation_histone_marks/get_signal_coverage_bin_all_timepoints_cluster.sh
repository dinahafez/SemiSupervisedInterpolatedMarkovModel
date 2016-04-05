#!/bin/bash
#$ -l h_vmem=20G
#$ -m bes
#$ -cwd
#$ -e error_jamm.txt
#$ -o error_jamm.txt


qsub run_BedToPeakCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/all_pos_14_stage_2_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/all_pos_14_stage_2_mid.bed /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/all_pos_14_stage_2 160 10




qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/all_pos_14_stage_2_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/all_pos_14_stage_2_mid.bed /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/all_pos_14_stage_2 160 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/all_pos_14_stage_2_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me1/H3K4me1.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/all_pos_14_stage_2_mid.bed /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/all_pos_14_stage_2 160 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/all_pos_14_stage_2_intersect_histone.reads /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K27ac/H3K27ac.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/all_pos_14_stage_2_mid.bed /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/all_pos_14_stage_2 160 10

#Background
qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/DHS_no_overlap_TSS_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_no_overlap_TSS_not_pos_14_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/DHS_no_overlap_TSS_not_pos_14 160 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/DHS_no_overlap_TSS_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me1/H3K4me1.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_no_overlap_TSS_not_pos_14_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/DHS_no_overlap_TSS_not_pos_14 160 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/DHS_no_overlap_TSS_intersect_histone.reads /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K27ac/H3K27ac.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_no_overlap_TSS_not_pos_14_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/DHS_no_overlap_TSS_not_pos_14 160 10




#Promoter
qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/DHS_overlap_TSS_intersect_histone.reads   /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_overlap_TSS_2bp_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me3/DHS_overlap_TSS 160 10 

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/DHS_overlap_TSS_intersect_histone.reads   /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me1/H3K4me1.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_overlap_TSS_2bp_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K4me1/DHS_overlap_TSS 160 10 

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/DHS_overlap_TSS_intersect_histone.reads  /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K27ac/H3K27ac.4.12.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/DHS_overlap_TSS_2bp_mid.bed  /data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Validation_Histone_marks/H3K27ac/DHS_overlap_TSS 160 10 




#################

qsub run_BedToPeakSignalCoverageMata_strandSpecific.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K27ac/DHS_overlap_TSS_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K27ac/H3K27ac.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/DHS_overlap_TSS_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K27ac/DHS_overlap_TSS 167 10 /data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed

qsub run_BedToPeakSignalCoverageMata_strandSpecific.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me1/DHS_overlap_TSS_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me1/H3K4me1.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/DHS_overlap_TSS_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me1/DHS_overlap_TSS 162 10 /data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed

qsub run_BedToPeakSignalCoverageMata_strandSpecific.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me3/DHS_overlap_TSS_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/DHS_overlap_TSS_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me3/DHS_overlap_TSS 160 10 /data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed


#background
qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me1/background_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me1/H3K4me1.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/background_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me1/background 162 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me3/background_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K4me3/H3K4me3.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/background_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K4me3/background 160 10

qsub run_BedToPeakSignalCoverageMata.sh /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K27ac/background_5kb_reads.bed /data/ohler/Dina/Drosophila/Research/data/Histone_modification_chip_seq/Jamm_peaks/H3K27ac/H3K27ac.filtered.peaks.NarrowPeaks /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/background_mid.bed /data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/Validation_Histone_marks/H3K27ac/background 167 10