library("ROCR")
args <- commandArgs(trailingOnly = TRUE)
inputFile="/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/5_mer/16_pos_single_ubiquitous_single.predictions";
  #args[1];
outputFile="";args[2];
inputFile;
outputFile;
currPreds <- read.table(inputFile, header = TRUE);
attach(currPreds);
pred <- prediction(currPreds$predictions,currPreds$labels);
perf <- performance(pred,"tpr","fpr")
plot(perf,col="black",lty=3, lwd=3)
