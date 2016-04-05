average_1=NULL;
average_2=NULL;

setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Baseline/")
sum = 0;
total=0;
allAUC=NULL;
for(i in 11:39)
{
  file3 = paste(i,"_DHS_closest_single_ubiquitous_single_Classifier.ROC", sep="")
 # print (file3)
  if (file.exists(file3))
  { 
    auc = read.table(file3, header=FALSE)$V3
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
  
  file3 = paste(i,"_pos_single_ubiquitous_single_Classifier.ROC", sep="")
  
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC2 = c(allAUC2,auc2)
       
  }
}

####################################################
#setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection//ModelStage1/MC_multi_output//ModelStage2/MC_multi_output//MergedDHS/5_mer/")

allAUC3=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file3 = paste(i,"_50kb_single_ubiquitous_single_Classifier.ROC", sep="")
  
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC3 = c(allAUC3,auc2)
   
    
  }
}



#########PLOT
#pdf("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/Graphs/Baselines.pdf")

par(ann=FALSE)
plot(jitter(rep(0,length(allAUC2)),amount=0.25), allAUC2,
     xlim=range(-0.5,2.5), ylim=range(0.4,1),
     axes=FALSE,frame.plot=TRUE, lwd=2)
points(jitter(rep(1,length(allAUC)), amount=0.25), allAUC, col=2, lwd=2)
points(jitter(rep(2,length(allAUC3)), amount=0.25), allAUC3, col=4, lwd=2)


##Add in the y-axis
axis(2, seq(0,1,by=0.1))

##Add in the x-axis labels
mtext("closestGene", side = 1, at=-0,line=1, cex = 1.3)
mtext("oneClosestDHS", side = 1, at=1,line=1, cex = 1.3)
#mtext("stage2_12", side = 1, at=2)
mtext("allDHS50kb", side = 1, at=2.0,line=1, cex = 1.3)
mtext("auROC", side=2, line = 3, cex = 1.3)
##Add in the means
segments(-0.25, mean(allAUC2), 0.25, mean(allAUC2))
segments(0.75, mean(allAUC), 1.25, mean(allAUC))
segments(1.75, mean(allAUC3), 2.25, mean(allAUC3))
#segments(1.75, mean(allAUC4), 2.25, mean(allAUC4))
#segments(2.75, mean(allAUC5), 3.25, mean(allAUC5))
#title("l1-logreg classification against ubiquitous")
#dev.off()
