#run scotts scripts

for i in {11..39..1}
do  

	file="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/"$i"_initialize_redfly_uniq_single.bed";
	output="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$i"_initialize_redfly_mid.bed";
	perl /data/ohler/Dina/Drosophila/Research/code/validaiton_histone_marks/DHScenterWindowToStdout.pl $file 2000 2000 > $output;
done

perl /data/ohler/Scott/Perl_Scripts/DHScenterWindowToStdout.pl "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/ubiquitous.bed" 2000 2000 > "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/ubiquitous_mid.bed";

for i in {11..39..1}
do  

	file="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/MC_5order/MC_ubiq_output/"$i"_pos.bed";
	output="/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/"$i"_pos_mid.bed";
	perl /data/ohler/Scott/Perl_Scripts/DHScenterWindowToStdout.pl $file 2000 2000 > $output;
done

perl /data/ohler/Scott/Perl_Scripts/DHScenterWindowToStdout.pl "/data/ohler/Dina/Drosophila/Research/data/DHS/DHS_overlap_TSS_10bp_withStrand.bed" 2000 2000 > "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Validation_Histone_marks/DHS_overlap_TSS_10bp_mid.bed";
