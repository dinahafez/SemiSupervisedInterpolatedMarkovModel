#read bed files 
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
permute = args[1];
enhancers_bed <- file.path("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/", "partnerTogene_intersect_DHS_sig_uniq.bed")
file1 = paste ("DHStoRandomeGene", permute, sep="_");
#file2 = paste(file1, ".bed_intersect_sexton_partnerTogene.bed.filtered.uniq.bed", sep="")
file2 = "model_intersect_sexton_partnerTogene_sig_uniq.bed"
mc_bed <- file.path("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene", file2)

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

file5 = paste ("Random_overlaps_hic_genes", permute, sep="_")
file5 = "model_overlaps_HiC_genes.bed";
file6 = paste ("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/", file5, sep="") 

file3 = paste ("Random_overlaps_hic", permute, sep="_")
file3 = "model_overlaps_HiC.bed";
file4= paste ("/data/ohler/Dina/Drosophila/Research/results_idr_0.2/AssignDHSToRandomGene/", file3, sep="")

for (i in 1:length(l1_uniq))
  {
  ind = sapply(l1,function(x, name){name %in% x}, l1_uniq[i])
    ind2 = sapply(l2,function(x, name){name %in% x}, l1_uniq[i])

    overlap = findOverlaps( mc_ranged[ind2,], enhancer_ranged[ind,], type='any',select='all')
    sub = subsetByOverlaps(  mc_ranged[ind2,],enhancer_ranged[ind,], type='any')
  
   
    countEnhancers = countEnhancers + nrow(sub);
 
    if (nrow(sub) > 0)
    {
      
      write.table(l1_uniq[i], file = file6,  
                  append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) #overlap_stark_gene_names_common_5.bed",
      
      Genes = c(Genes, drop(sub[,'name'])[[1]]);
     
      total= nrow(sub);
      idx <- unique(subjectHits(overlap))
      for (k in 1:total)
      {
        
        aline= c( toString(sub[['space']][k]), start(sub[['ranges']][k]), end(sub[['ranges']][k]),  l1_uniq[i]);
       
        FF <- as.matrix(t(aline))
    
        write.table(FF, file =file4 ,
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) #overlap_stark_common_5.bed", 
      }
    }
    
    
}
Genes_uniq = unique(Genes)
