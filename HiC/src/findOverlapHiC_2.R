#read bed files 
library(rtracklayer)
#dir = "/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap"
enhancers_bed <- file.path("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/", "partnerTogene_intersect_model_uniq.bed")
#mc_bed <- file.path("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output", "all_clusters_stage_2.bed")

for (i in c(14,12,10) )
{
	for (j in c(14,10)) 
	{
		file1 = paste("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2_pos_", i, sep="");
		file2 = paste(file1,"all_pos_", sep="/")
		file3 = paste(file2,j,sep="")
		file4 = paste(file3,"_uniq_genes.bed", sep="")
		print (file4)
mc_bed <- file.path(file4);#"all_pos_5_for_stark_overlap.bed")

enhancer <- import(enhancers_bed)
mc <- import(mc_bed)

enhancer_ranged = as(enhancer, "RangedData")
mc_ranged = as(mc, "RangedData")

#get names of genes from enhancers
l1 = drop( enhancer_ranged[,'name'])[[1]]
l2 = drop( mc_ranged[,'name'])[[1]]
l1_uniq = unique(l1)


#should loop on all gene names in l1
countEnhancers=0;
Genes=NULL;

for (i in 1:length(l1_uniq))
  {
  ind = sapply(l1,function(x, name){name %in% x}, l1_uniq[i])
    ind2 = sapply(l2,function(x, name){name %in% x}, l1_uniq[i])

    overlap = findOverlaps( mc_ranged[ind2,], enhancer_ranged[ind,], type='any',select='all')
  #  sub = subsetByOverlaps(  enhancer_ranged[ind,],mc_ranged[ind2,], type='any')
    
    sub = subsetByOverlaps(  mc_ranged[ind2,],enhancer_ranged[ind,], type='any')
  
   
    countEnhancers = countEnhancers + nrow(sub);
    if (nrow(sub) > 0)
    {
	fileGenes = paste(file3,"_overlap_HiC_genes.bed", sep="")
	print (fileGenes);
      write.table(l1_uniq[i], file = fileGenes,
                  append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) #overlap_stark_gene_names_common_5.bed",
      
      Genes = c(Genes, drop(sub[,'name'])[[1]]);
      #countGenes = countGenes + length(overlap);
      total= nrow(sub);
      idx <- unique(subjectHits(overlap))
      for (k in 1:total)
      {
        
        #aline= c( names(sub)[k], start(sub[['ranges']][k]), end(sub[['ranges']][k]),  l1_uniq[i]);
        aline= c( toString(sub[['space']][k]), start(sub[['ranges']][k]), end(sub[['ranges']][k]),  l1_uniq[i]);
       
        FF <- as.matrix(t(aline))
       file_out = paste(file3,"_overlap_HiC.bed", sep="" )
        write.table(FF, file = file_out ,
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) #overlap_stark_common_5.bed", 
      }
    }
    
    
}
Genes_uniq = unique(Genes)
}}
