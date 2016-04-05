average_1=NULL;
average_2=NULL;
####################################################
####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//5_mer/")
setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/MergedDHS/5_mer/")
sum = 0;
total=0;
allAUC=NULL;
for(i in 15:39)
{

  file = paste(i,"_pos_10_single_merged_ubiquitous_all_single_merged_Classifier.ROC", sep="")
  
  if (file.exists(file))
  { 
    auc = read.table(file, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC = c(allAUC,auc2)
   
    print (file)
  
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
allAUC2=NULL;
sum = 0;
total=0;
for(i in 15:39)
{
 
  file = paste(i,"_pos_10_single_merged_stage2_all_ubiquitous_all_single_merged_stage2_all_Classifier.ROC", sep="")
  
  if (file.exists(file))
  {
    auc = read.table(file, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC2 = c(allAUC2,auc2)
    print (file)
  }
}

average = sum/total;
average_1=c(average_1,average);

###################
allAUC3=NULL;
sum = 0;
total=0;
for(i in 15:39)
{
  
  file = paste(i,"_TSS_all_single_merged_ubiquitous_TSS_all_single_merged_Classifier.ROC", sep="")
  
  if (file.exists(file))
  {
    auc = read.table(file, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC3 = c(allAUC3,auc2)
    print (file)
  }
}


#########PLOT
#pdf("/Users/Dina/Google Drive/Thesis_v1.1/Pictures/mc_roc.pdf")

par(ann=FALSE)
plot(jitter(rep(0,length(allAUC)),amount=0.25), allAUC,
     xlim=range(-0.5,2.5), ylim=range(0.5,1),
     axes=FALSE,frame.plot=TRUE, lwd=2, col='blue')
points(jitter(rep(1,length(allAUC3)), amount=0.25), allAUC3, col='black', lwd=2)
points(jitter(rep(2,length(allAUC2)), amount=0.25), allAUC2, col="red", lwd=2)
#points(jitter(rep(3,length(allAUC4)), amount=0.25), allAUC4, col=4)

##Add in the y-axis
axis(2, seq(0,1,by=0.1), cex=1.5)


##Add in the x-axis labels
mtext("Distal-only", side = 1, at=0, cex = 1.5, line=1)
mtext("Distal+TSS", side = 1, at=2, cex = 1.5, line=1)
mtext("TSS-only", side = 1, at=1, cex = 1.5, line=1)

mtext("auROC", side=2, line = 3,  cex = 1.5)
##Add in the means
segments(-0.25, mean(allAUC), 0.25, mean(allAUC))
segments(0.75, mean(allAUC3), 1.25, mean(allAUC3))
segments(1.75, mean(allAUC2), 2.25, mean(allAUC2))
#segments(2.75, mean(allAUC4), 3.25, mean(allAUC4))
title(" ")
#dev.off()
