setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/")

figfile = paste("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/Graphs/histograms/overlap_","cdf_2.pdf", sep="");
#pdf(figfile)
allOver=NULL
#par(mfrow=c(3,4)) 
for(i in 15:38)
{
#  file = paste("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/",i, sep="")
  file2 = paste(i,"_overlap.bed",sep="")
  print (file2)
  
  if (file.exists(file2))
  {
      overlap = read.table(file2, h=F)
      j=i-10
      t = paste("cluster ",j)
      t2 = paste(t,"R",sep="")
      allOver = c(allOver,overlap$V5)
      #hist(overlap$V5, breaks = 100, main=t2, xlab="#common DHSs")
      #plot(ecdf(overlap$V5),  xlab='#common DHSs',ylab="CDF", main=t2)
     # t = paste("cluster ",i)
    #  title(t)
     
  }
  
}
pdf("~/Google Drive/Thesis_v1.1/Pictures/overlap_all.pdf")
hist(allOver, breaks = 50, main="", xlab="shared DHSs", col="grey", cex.lab=1.5)
dev.off()
#dev.off()