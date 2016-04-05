#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error.out -j y
#$ -l h_vmem=20G
#$ -l mem_free=20G
#$ -e error.log

########  understand data

#Pre: all files for initialization for each cluster have been generated, labeled, unlabeled, and background files, and saved in result directory
dataDrive="/Users/Dina/Duke/Uwe_Lab/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/";


#############################   semi-supervised learning -- MC with EM    ###############################

MC_order=$1;	 #Mc order
clusterNo=$2	 #cluster no to process
kmer=$3;         #no of kmers
background="ubiq";   #could be ubiq or other

				
#create directories 
MC_order_path="MC_"$MC_order"order/";  
parentDir=$dataDrive;
fullPath=$parentDir$MC_order_path;
if [ ! -d "$fullPath" ]
then
	mkdir $fullPath;
fi

if [ $background == "other" ]
then
	
	fullPath=$fullPath"MC_other_output/";
else
	fullPath=$fullPath"MC_ubiq_output/";
fi

if [ ! -d "$fullPath" ]
then
	mkdir  $fullPath;
fi
if [ ! -d "$fullPath$kmer"_mer"" ]
then
	mkdir  $fullPath$kmer"_mer";
fi
if [ ! -d "$fullPath$kmer"_mer/tempProcessingFiles"" ]
then
	mkdir  $fullPath$kmer"_mer/tempProcessingFiles";
fi

#run Makov Chain 
dhsLabeledFileName=$dataDrive$clusterNo"_initialize_single.fa"; 
dhsUnLabeledFileName=$dataDrive$clusterNo"_unlabeled_single.fa";
#Markov chain background could be other or ubiq
if [ $background == "other" ]
then
	dhsBackgroundFileName=$parentDir$clusterNo"_neg.fa";
	perl getNegativeFiles.pl $parentDir  $clusterNo;
else
	dhsBackgroundFileName=$parentDir"ubiquitous_single.fa";  
fi
posFile=$clusterNo"_pos_single";
dhsLabeledPositiveFileName=$fullPath$posFile".fa";

if [ ! -s $dhsLabeledPositiveFileName ]
then
	java -Xmx7g -jar /Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/EMwithMC/target/EMwithMC-1.4-SNAPSHOT-jar-with-dependencies.jar $dhsLabeledFileName $dhsUnLabeledFileName $dhsBackgroundFileName $dhsLabeledPositiveFileName $MC_order "3";

fi

############################     get kmer counts    #################################

echo "Generating kmers ";

if [ ! -s $fullPath$kmer"_mer/"$posFile".tab" ]
then 
	python /Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/sparseClassifier/getKmerCount/src/getKmerCount.py $fullPath $posFile $fullPath $kmer
fi

ubiqFile="ubiquitous_single";
if [ ! -s $fullPath$kmer"_mer/"$ubiqFile".tab" ]
then
	python /Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/sparseClassifier/getKmerCount/src/getKmerCount.py  $fullPath $ubiqFile $fullPath $kmer
fi

echo "Generating kmers finished";

############################sparse classification   ##################################
#classify between positive vs ubiquitous

echo "running classifier";

ClassifyFile=$fullPath$kmer"_mer/"$posFile"_"$ubiqFile"_Classifier.ROC";
echo $ClassifyFile;
if [ ! -s $ClassifyFile ]
then	
	java -Xmx7g -jar /Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/ClassifyKmerCounts/target/ClassifyKmerCounts-1.4-SNAPSHOT-jar-with-dependencies.jar $fullPath $posFile $ubiqFile $kmer
fi

###################     plot ROC curve #######################

Rscript /Users/Dina/Duke/Uwe_Lab/Drosophila/Research/code/R_plots/plot_roc.R $fullPath$kmer"_mer/"$posFile"_ubiquitous_single.predictions" $fullPath$kmer"_mer/"$posFile"_ubiquitous_single.pdf" 



