#! /bin/bash

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
Rscript call_adapteds.R -u --exclude=CLM\|FLC4\|Stan --gens=8 --cutoff=0.05 --base_name=cutoff-5 --outdir=../data/fitness_data/fitness_calls ../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv Ancestor_YPD_2N
