#read bed files 
library(rtracklayer)
enhancers_bed <- file.path("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/", "partnerTogene_intersect_DHS_genes_in_study_stage1.bed")
mc_bed <- file.path("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results", "stage1_14_intersect_partnerTogene_stage1.bed.filtered.uniq.bed")

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
      
      write.table(FF, file ="/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/HiC/results/stage1_14_overlap_hic.bed" ,
                  append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) 
    }
  }
  
  
}
