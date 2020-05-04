#!/bin/bash

# analysis pipeline for plt

###########################
# identify neutral BC set #
###########################

# call_adapteds.R
#
# inputs:
  # BFA fitness files
  # putatively neutral set identifier
  # stopping criteria for outlier detection
  # outdir name

# dBFA2:
Rscript call_adapteds.R -u \
  --exclude=CLM\|FLC4\|Stan \
  --gens=8 \
  --cutoff=0.05 \
  --base_name=dBFA2_cutoff-5 \
  --outdir=../data/fitness_data/fitness_calls \
  ../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv \
  Ancestor_YPD_2N

# hBFA1
Rscript call_adapteds.R -u \
  --exclude=X48Hr \
  --gens=8 \
  --cutoff=0.05 \
  --base_name=hBFA1_cutoff-5 \
  --outdir=../data/fitness_data/fitness_calls \
  ../data/fitness_data/fitness_estimation/hBFA1_s_03_23_18_GC_cutoff_5.csv \
  YPD_alpha

# hBFA2
# This assay is special because we do not have a neutral reference for 2N clones
# Rscript call_adapteds.R -u \
#   --exclude=CLM\|X08M_NaCl\|X48Hr \
#   --gens=8 \
#   --cutoff=0.05 \
#   --base_name=hBFA2_cutoff-5 \
#   --outdir=../data/fitness_data/fitness_calls \
#   ../data/fitness_data/fitness_estimation/hBFA2_s_03_23_18_GC_cutoff_5.csv \
#   YPD_alpha

# run autodiploid flags on hBFA1 adapted calls
Rscript filter_autodiploids.R \
  --use_iva \
  --exclude=X48 \
  --gens=8 \
  --cutoff=0 \
  --base_name=hBFA1_cutoff-5 \
  --outdir="../data/fitness_data/fitness_calls" \
  "../data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_2020-03-03.csv" \
  "autodiploids"


  if (debug_status == TRUE) {
      arguments <- list()
      arguments$use_iva     <- TRUE
      arguments$infile      <- "../data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_2020-03-03.csv"
      arguments$outdir      <- "../data/fitness_data/fitness_calls"
      arguments$autodip_tag <- 'autodiploids'
      arguments$base_name   <- 'hBFA1_cutoff-5'
      arguments$exclude     <- 'X48'
      arguments$cutoff      <- 0
      arguments$gens        <- 8



# run t-SNE cluster on all lineages (hBFA1 and dBFA2)
Rscript cluster_lineages.R -u \
  --gens=8 \
  --exclude=CLM\|FLC4\|Stan\|X48Hr \
  --outdir=../data/fitness_data/fitness_calls/clusters \
  ../data/fitness_data/fitness_calls/dBFA2_cutoff-5_adapteds_2020-03-03.csv \
  ../data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_autodips.csv
