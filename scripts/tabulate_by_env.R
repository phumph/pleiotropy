#! /usr/bin/env Rscript

# tabulate_by_env.R

# This script takes input fitness data (post-clustering)
# as well as mutation data linked with barcodes and
# generates counts of lineages across environments.
#
# Outputs comprise two Supplemental Tables.

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))

source(file.path('scripts/src/pleiotropy_functions.R'))



# ==== #
# main #
# ==== #

"tabulate_by_env.R

Usage:
    cluster_lineages.R [--help]
    cluster_lineages.R [options] <infile>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations

Arguments:
    infile                        Input file(s) containing fitness calls for BFA run.
" -> doc

# define default args for debug_status == TRUE
args <- list(
  use_iva = TRUE,
  infile = "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_autodips.csv",
  outdir = "data/fitness_data/fitness_calls",
  gens = "8",
  exclude = "X48Hr"
)

debug_status <- FALSE

cat("\n*********************\n")
cat("* tabulate_by_env.R *\n")
cat("*********************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")