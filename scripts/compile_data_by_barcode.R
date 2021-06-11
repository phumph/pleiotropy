#! /usr/bin/env Rscript
# compile_data_by_barcode.R

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggplot2)))

source(file.path("scripts/src/pleiotropy_functions.R"))


main <- function(arguments) {
  
  fitness_data <- read.table(arguments$fitness_file,
                             sep = ",",
                             header = T,
                             stringsAsFactors = F)
  
  cluster_data <- read.table(arguments$cluster_file,
                             sep = ",",
                             header = T,
                             stringsAsFactors = F)
  
  mutations_data <- read.table(arguments$mutations_file,
                               sep = ",",
                               header = T,
                               stringsAsFactors = F)
  
  fitness_data$has_wgs <- fitness_data$Diverse.BC %in% mutations_data$Diverse.BC
  
  fitness_data %>%
    dplyr::filter(!Subpool.Environment %in% c("none", "not_read")) %>%
    dplyr::left_join(dplyr::select(cluster_data, Full.BC, cluster),
                     by = "Full.BC") ->
    fitness_data
  
  fit_dat_cols <- c("is_adapted", "neutral_set", "has_wgs", "cluster")
  if (grepl("^hBFA", basename(arguments$fitness_file))) {
    fit_dat_cols <- c(fit_dat_cols, "autodip")
  }
  
  fitness_data %>%
    prep_fit_matrix(means_only = FALSE,
                    excludes = arguments$exclude,
                    iva_s = arguments$use_iva,
                    gens = arguments$gens) ->
    fitness_data_prepped
  
  fitness_data_prepped[[1]] %>%
    as.data.frame() %>%
    dplyr::mutate(Full.BC = row.names(.)) %>%
    tidyr::pivot_longer(cols = grep("_s$", names(.), value = T),
                        names_to = "bfa_env", values_to = "s") %>%
    dplyr::left_join(
      fitness_data[, c("Full.BC", "Subpool.Environment", fit_dat_cols)],
      by = "Full.BC") -> fitness_matr
  
  fitness_matr %>%
    normalize_envs() ->
    fitness_matr_norm
  
  # exclude autodiploids if hBFA
  if (grepl("hBFA", basename(arguments$fitness_file))) {
    fitness_matr_norm %>%
      dplyr::filter(autodip == FALSE) %>%
      dplyr::select(-autodip) ->
      fitness_matr_norm
  }
  
  # add stderr data
  fitness_data_prepped[[2]] %>%
    as.data.frame() %>%
    dplyr::mutate(Full.BC = row.names(.)) %>%
    tidyr::pivot_longer(cols = grep("_s_err", names(.), value = T),
                        names_to = "bfa_env", values_to = "s_se") %>%
    dplyr::left_join(
      fitness_data[, c("Full.BC", "Subpool.Environment", fit_dat_cols)],
      by = "Full.BC") %>%
    normalize_envs() %>%
    dplyr::select(Full.BC, bfa_env, s_se, source)->
    fitness_matr_norm_err
  
  fitness_matr_norm %>%
    dplyr::left_join(fitness_matr_norm_err, by = c("Full.BC", "bfa_env", "source")) %>%
    dplyr::select(Full.BC, source, cluster, bfa_env, s, s_se, is_adapted, neutral_set, has_wgs) ->
    fitness_matr_norm_comb
  
  bfa_prfx <- strsplit(basename(arguments$fitness_file), "_adapteds")[[1]][1]
  
  fitness_matr_norm_comb %>%
    write_out(out_dir = file.path(arguments$outdir),
              base_name = basename(arguments$fitness_file),
              str_to_append = "_compiled_data_by_barcode")
  
}

# ==== #
# main #
# ==== #

"compile_data_by_barcode.R

Usage:
    compile_data_by_barcode.R [--help]
    compile_data_by_barcode.R [options] <fitness_file> <cluster_file> <mutations_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./data/compiled]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations

Arguments:
    fitness_file                  Input from adapted calling routine
    cluster_file                  Input with cluster IDs (from cluster_lineages.R)
    mutations_file                Input with mutation data (from combine_BCs_and_WGS.R)
" -> doc

# define default args for debug_status == TRUE
args <- list(
  use_iva = TRUE,
  fitness_file = "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_autodips.csv",
  cluster_file = "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapted_w_clusts.csv",
  mutations_file = "data/mutation_data/mutations_by_bc.csv",
  outdir = "data/combined",
  gens = "8",
  exclude = "X48Hr"
)

debug_status <- FALSE

cat("\n***************************\n")
cat("* compile_data_by_barcode *\n")
cat("***************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
