#! /usr/bin/env Rscript

# tabulate_by_env.R

# This script takes input fitness data (post-clustering)
# as well as mutation data linked with barcodes and
# generates counts of lineages across environments.
#
# Outputs comprise two Supplemental Tables.

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggplot2)))

source(file.path('scripts/src/pleiotropy_functions.R'))

# ------------- #
# function defs #
# ------------- #

normalize_envs <- function(df,
                           remove_str_from_source = "_alpha$",
                           remove_str_from_bfa_env = "^X") {
  
  stopifnot(is.data.frame(df) & all(c("source", "bfa_env") %in% names(df)))
  
  df$source <- sapply(df$source, function(x) gsub(remove_str_from_source, "", x))
  df$bfa_env <- sapply(df$bfa_env, function(x) gsub(remove_str_from_bfa_env, "", x))
  
  return(df)
}


flag_home <- function(df) {
  
  df$home <- df$bfa_env == df$source
  
  return(df)
}


flag_pleiotropy <- function(df, z_cut = 1.96,
                            source_col = "source",
                            cluster_col = "cluster",
                            home_col = "home") {
  
  df$fitness_diff <- ((abs(df$s / df$s_se)) > 1.96)
  df$pleio_pos <- df$fitness_diff == TRUE & df$s > 0
  df$pleio_neg <- df$fitness_diff == TRUE & df$s < 0
  
  df %>%
    dplyr::filter(!!sym(home_col) == FALSE) %>%
    dplyr::group_by(!!sym(source_col),
                    !!sym(cluster_col)) %>%
    dplyr::summarise(s_mu_away = mean(s),
                     n_bfa_envs = length(unique(bfa_env)),
                     n_pleio_pos = sum(pleio_pos),
                     n_pleio_neg = sum(pleio_neg),
                     n_pleio_tot = sum(fitness_diff)) %>%
    dplyr::left_join((df %>%
                        dplyr::filter(!!sym(home_col) == TRUE) %>%
                        dplyr::group_by(!!sym(source_col),
                                        !!sym(cluster_col)) %>%
                        dplyr::summarise(s_mu_home = mean(s))),
                     by = c(source_col, cluster_col)) ->
    df_summarized
    
  return(df_summarized)
}


plot_pleiotropy <- function(dat) {
  
  dat %>%
    ggplot(aes(x = s_mu_home, y = n_pleio_neg)) +
    geom_point() +
    facet_wrap(~source)
  
  dat %>%
    ggplot(aes(x = s_mu_home, y = s_mu_away))
}



main <- function(arguments) {
  
  infile <- read.table(arguments$infile,
                       sep = ",",
                       header = T,
                       stringsAsFactors = F)
  
  infile %>%
    normalize_envs() %>%
    flag_home() %>%
    flag_pleiotropy ->
    dat
  
  
}

# ==== #
# main #
# ==== #

"summarise_clusters.R

Usage:
    summarise_clusters.R [--help]
    summarise_clusters.R [options] <infile>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]

Arguments:
    infile                        Input file(s) containing cluster means from cluster_lineages.R.
" -> doc

# define default args for debug_status == TRUE
args <- list(
  infile = "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapted_w_clust_means.csv",
  outdir = "data/fitness_data/fitness_calls"
)

debug_status <- TRUE

cat("\n************************\n")
cat("* summarise_clusters.R *\n")
cat("************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
