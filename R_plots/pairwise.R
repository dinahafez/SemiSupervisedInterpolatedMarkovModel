average_1=NULL;
average_2=NULL;

setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/MergedDHS/pairwise/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file3 = paste("15_pos_10_single_merged_",i,sep="")
  file4 = paste(file3,"_pos_10_single_merged_Classifier.ROC", sep="")
 # print (file3)
  if (file.exists(file4))
  { 
    auc = read.table(file4, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC = c(allAUC,auc2)

    print (file3)
    
  }
}


####################################################
#setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output//ModelStage2/MC_multi_output//MergedDHS/5_mer/")

allAUC2=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file3 = paste("16_pos_10_single_merged_",i,sep="")
  file4 = paste(file3,"_pos_10_single_merged_Classifier.ROC", sep="")
  # print (file3)
  if (file.exists(file4))
  { 
    auc = read.table(file4, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC2 = c(allAUC2,auc2)
       
  }
}

average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output//ModelStage2/MC_multi_output//MergedDHS/5_mer/")

allAUC3=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file3 = paste(i,"_pos_12_single_merged_ubiquitous_all_single_merged_Classifier.ROC", sep="")
  
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC3 = c(allAUC3,auc2)
   
    
  }
}

####################################################
#setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output//MergedDHS/5_mer/")

allAUC4=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file3 = paste("17_pos_10_single_merged_",i,sep="")
  file4 = paste(file3,"_pos_10_single_merged_Classifier.ROC", sep="")
  # print (file3)
  if (file.exists(file4))
  { 
    auc = read.table(file4, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC4 = c(allAUC4,auc2)
  
  }
}


####################################################
#setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output//ModelStage2/MC_multi_output//MergedDHS/5_mer/")

allAUC5=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file3 = paste("33_pos_10_single_merged_",i,sep="")
  file4 = paste(file3,"_pos_10_single_merged_Classifier.ROC", sep="")
  # print (file3)
  if (file.exists(file4))
  { 
    auc = read.table(file4, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC5 = c(allAUC5,auc2)
    #print (file3);
  }
}


#########PLOT
#pdf("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/Graphs/overlap_stage_2_cutoffs_mergedDHS.pdf")

par(ann=FALSE)
plot(jitter(rep(0,length(allAUC)),amount=0.25), allAUC,
     xlim=range(-0.5,1.5), ylim=range(0.45,1),
     axes=FALSE,frame.plot=TRUE, lwd=2)
points(jitter(rep(1,length(allAUC2)), amount=0.25), allAUC2, col=3, lwd=2)
#points(jitter(rep(2,length(allAUC3)), amount=0.25), allAUC3, col=4, lwd=2)
#points(jitter(rep(2,length(allAUC4)), amount=0.25), allAUC4, col=4, lwd=2)
#points(jitter(rep(3,length(allAUC5)), amount=0.25), allAUC5, col=2, lwd=2)

##Add in the y-axis
axis(2, seq(0,1,by=0.1))

##Add in the x-axis labels
mtext("cluster_18R", side = 1, at=-0,line=1)
mtext("cluster_19R", side = 1, at=1,line=1)
#mtext("stage2_12", side = 1, at=2)
#mtext("cluster_7R", side = 1, at=2.0,line=1)
#mtext("cluster_23R", side = 1, at=3.0, line=1)
mtext("auROC", side=2, line = 3)
##Add in the means
segments(-0.25, mean(allAUC), 0.25, mean(allAUC))
segments(0.75, mean(allAUC2), 1.25, mean(allAUC2))
#segments(1.75, mean(allAUC3), 2.25, mean(allAUC3))
segments(1.75, mean(allAUC4), 2.25, mean(allAUC4))
#segments(2.75, mean(allAUC5), 3.25, mean(allAUC5))
#title("l1-logreg classification against ubiquitous")
#dev.off()
