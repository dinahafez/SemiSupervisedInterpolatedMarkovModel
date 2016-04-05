average_1=NULL;
average_2=NULL;

Stark_dir="/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/5_mer/";
Stark_dir_stage2="/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output/5_mer/";
####################################################
####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//5_mer/")
setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.1/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/5_mer")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{

  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  
  if (file.exists(file3))
  { 
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC = c(allAUC,auc2)
    sum = sum+auc;
    total=total+1;
    print (file3)
  
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//MC_2nd_3order/MC_ubiq_output/5_mer/")
setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.1/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/MC_ubiq_output/5_mer")
allAUC2=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
 
    file = paste(i,"_pos_all_pos_not_", sep="")
    file2 = paste(file,i,sep="")
    file3 = paste(file2,"_Classifier.ROC",sep="")
    file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
       
  file_stark=paste(Stark_dir_stage2,file3,sep="")
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC2 = c(allAUC2,auc2)
    sum = sum+auc;
    total=total+1;
    print (file3)
  }
}

average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//5_mer/")
setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/5_mer")
allAUC3=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  
  file_stark=paste(Stark_dir,file3,sep="")
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC3 = c(allAUC3,auc2)
    sum = sum+auc;
    total=total+1;
  
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//MC_2nd_3order/MC_ubiq_output/5_mer/")
setwd("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/MC_ubiq_output/5_mer")
allAUC4=NULL;
sum = 0;
total=0;
for(i in 11:39)
{

  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  
  
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    auc2 = as.numeric(as.character(auc))
    allAUC4 = c(allAUC4,auc2)
    sum = sum+auc;
    total=total+1;
  }
}

average = sum/total;
average_1=c(average_1,average);

#########PLOT
pdf("/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/Graphs/Initialization_both_diff_idr_cutoffs_MC_order_2.pdf")

par(ann=FALSE)
plot(jitter(rep(0,length(allAUC)),amount=0.25), allAUC,
     xlim=range(-0.5,3.5), ylim=range(0.4,1),
     axes=FALSE,frame.plot=TRUE, lwd=2)
points(jitter(rep(1,length(allAUC2)), amount=0.25), allAUC2, col=4, lwd=2)
points(jitter(rep(2,length(allAUC3)), amount=0.25), allAUC3, col=3, lwd=2)
points(jitter(rep(3,length(allAUC4)), amount=0.25), allAUC4, col=2, lwd=2)

##Add in the y-axis
axis(2, seq(0,1,by=0.1))

##Add in the x-axis labels
mtext("idr_0.1_MC_3", side = 1, at=-0.25)
mtext("idr_0.1_MC_2", side = 1, at=1)
mtext("idr_0.2_MC_3", side = 1, at=2)
mtext("idr_0.2_MC_2", side = 1, at=3.0)
##Add in the means
segments(-0.25, mean(allAUC), 0.25, mean(allAUC))
segments(0.75, mean(allAUC2), 1.25, mean(allAUC2))
segments(1.75, mean(allAUC3), 2.25, mean(allAUC3))
segments(2.75, mean(allAUC4), 3.25, mean(allAUC4))
title("Init both, different IDR cutoffs & MC order")
dev.off()
