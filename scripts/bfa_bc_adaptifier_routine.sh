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

# next step: post-process hBFA1 adapted columns to flag auto-diploids, both adapted and ancestral.
