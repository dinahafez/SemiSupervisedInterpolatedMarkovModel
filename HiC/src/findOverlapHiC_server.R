#read bed files 
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
enhancers_bed <- file.path("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/", "partnerTogene_intersect_DHS_genes_in_study_stage1.bed")
mc_bed <- file.path("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/", args[1])


enhancer <- import(enhancers_bed)
mc <- import(mc_bed)

enhancer_ranged = as(enhancer, "RangedData")
mc_ranged = as(mc, "RangedData")

#get names of genes from enhancers
l1 = drop( enhancer_ranged[,'name'])[[1]]
l2 = drop( mc_ranged[,'name'])[[1]]
l1_uniq = unique(l2)


#should loop on all gene names in l1
countEnhancers=0;
Genes=NULL;

for (i in 1:length(l1_uniq))
{
  ind = sapply(l1,function(x, name){name %in% x}, l1_uniq[i])
  ind2 = sapply(l2,function(x, name){name %in% x}, l1_uniq[i])
    
  sub = subsetByOverlaps(  mc_ranged[ind2,],enhancer_ranged[ind,], type='any')
  
  if (nrow(sub) > 0)
  {
    total= nrow(sub);
    for (k in 1:total)
    {
      aline= c( toString(sub[['space']][k]), start(sub[['ranges']][k]), end(sub[['ranges']][k]),  l1_uniq[i]);
      
      FF <- as.matrix(t(aline))
      file1 = paste ("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/",args[1], sep="");   
	file2 = paste (file1, ".overlap.hic.bed", sep="");
       write.table(FF, file = file2,
                  append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) 
    }
  }
  
  
}
