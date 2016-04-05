file=$1
cd "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene"
perl ../../code/MC_AnaylzeOutput/convertFastaToBed.pl $file 
cut -f1,2,3,5 $file".bed"| sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq > $file"_gene.bed"
intersectBed -a $file"_gene.bed" -b "partnerTogene_intersect_DHS_genes_in_study_stage1.bed" -u | uniq > $file"_gene_intersect_partner_Togene.bed"

perl filter_genes_in_HiC.pl $file"_gene_intersect_partner_Togene.bed"
#Rscript findOverlapHiC_server.R $file"_gene_intersect_partner_Togene.bed.filtered.uniq.bed" 
