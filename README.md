# SemiSupervisedInterpolatedMarkovModel

This repository contains all the work done for the second project in m y PhD

The model links enhancers to their target genes based on their sequences.

DHSs are first alinged, peaks are called

AssignDHSToClusters_Jamm: Contains Perl scripts that assign DHSs to genes first based on old ways, then it initialized the assignment with known data

EMwithMC: Java code that implements the model

ClassifyKmerCounts: Java code. To validate the model output, we run a logisitc regression. This folder contains the code to do classification

HiC: All scirpts used to do the analysis with the HiC data

validation_TF_enrichment: Scripts for TF validation 

validation_histone_marks: Scripts for histone marks enrichment
