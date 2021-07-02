#!/usr/bin/env Rscript

# This script is adapted from the R scripts developed by Fernando Pérez Villatoro [https://github.com/farkkilab/ovaHRDscar]
# Implementation of NtAI, LST and HRD-LOH in R 

#########################################################################
############################     HOW TO USE    ##########################
#########################################################################

# 'x' is a matrix of segmented output from ASCAT :
# 1: sample_id
# 2: chromosome (numeric)
# 3: segment_start
# 4: segment_end
# 5: total_cn
# 6: nA
# 7: nB

# Each row represent a ASCN segment and total_cn = sum of nA and nB
# Note: For each segment, the package will re-order the A and B copy number values, considering as A the one with highest ASCN.
# get.ovaHRDscars("F:/Documents/scarHRD/examples/segments.txt",chrominfo = "grch37")

# Running parameters:
# seg -- Input file name
# reference -- the reference genome used, grch38 or grch37 (default: grch38)
library(data.table)

source("/home/vijay/my_projects/ovaHRDscar/R/calc.ai_new.R")
source("/home/vijay/my_projects/ovaHRDscar/R/features_LST.R")
source("/home/vijay/my_projects/ovaHRDscar/R/get_HRDs.R")
source("/home/vijay/my_projects/ovaHRDscar/R/get_LOH.R")
source("/home/vijay/my_projects/ovaHRDscar/R/preprocess.hrd.R")
source("/home/vijay/my_projects/ovaHRDscar/R/rm.chr.deletions.R")
source("/home/vijay/my_projects/ovaHRDscar/R/shrink.seg.ai.R")
source("/home/vijay/my_projects/ovaHRDscar/R/shrink.seg.ai.wrapper.R")

seg <- as.data.frame(read.table("/home/vijay/my_projects/HRD_pipeline/hrd_ascat_results/GSE118527/GSE118527_tumor_only/GSE118527_tumor_only_segments2.txt",header = T))
test <- get.ovaHRDscars(seg,chrominfo = "grch37")
test
write.table(test,"/home/vijay/my_projects/HRD_pipeline/hrd_ascat_results/GSE118527/GSE118527_tumor_only/GSE118527_HRDscores.txt",sep = '\t',col.names = T, quote = F)
