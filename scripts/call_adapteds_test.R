#! /usr/bin/env Rscript

# call_adapteds.R
#
# this script will take fitness files as input
# and handle a few parameters
# before outputting a table of barcodes and adapted status along with their fitness values.

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(docopt)))

# ------------------------- #
# setup command line inputs #
# ------------------------- #

'call_adapteds.R

Usage:
    call_adapteds.R [--help | --version]
    call_adapteds.R [options] <infile> <neutral_col>

Options:
    -h --help                     Show this screen.
    -v --version                  Show version.
    -o --outdir=<outdir>          Output directory [default: ./]
    -b --base_name=<name>         Base name for output file [default: tmp]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -c --cutoff=<pval>            P-value cutoff for outlier trimming [default: 0.001]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations
    -r --reps_iter=<sim_reps>     Number of replicate simulations for determining neutral set [default: 500]
    -f --reps_final=<final_reps>  Number of replicate simulations for final neutral probability distribution [default: 1000]

Arguments:
    infile                        Input file containing fitness calls for BFA run.
    neutral_col                   String denoting set of putatively neutral BCs
' -> doc

# -------------------- #
# function definitions #
# -------------------- #

run_args_parse <- function(debug_status) {

  if (debug_status == TRUE) {
    arguments <- list()
    arguments$use_iva     <- TRUE
    arguments$infile      <- "../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv"
    arguments$outdir      <- "../data/fitness_data/fitness_calls"
    arguments$neutral_col <- 'Ancestor_YPD_2N'
    arguments$name        <- 'cutoff-5'
    arguments$reps_iter   <- 500
    arguments$reps_final  <- 1000
    arguments$exclude     <- 'CLM|FLC4|Stan'
    arguments$cutoff      <- 0.05
    arguments$gens        <- 8
    arguments$base_name   <- 'dBFA2'
  } else if (debug_status == FALSE) {
    arguments <- docopt(doc, version = 'call_adapteds v.1.0')
  }
  return(arguments)
}

# ---- #
# MAIN #
# ---- #

# run
debug_status <- FALSE
arguments <- run_args_parse(debug_status)
print(arguments)
