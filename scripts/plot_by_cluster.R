#! /usr/bin/env Rscript

# plot_by_cluster.R

# This script takes input from summarise_clusters.R
# and plots fitness traces and pleiotropy patterns
# accross environments

# Outputs 
# 1. individual plots per cluster per source environment
# 2. summary plot of pleiotropic effects per source environment

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))

source(file.path("scripts/src/pleiotropy_functions.R"))

# ------------- #
# function defs #
# ------------- #

rescale_fitness <- function(df, bfa_envs, scale_factor, fit_col = "s") {
  df[[fit_col]][df$bfa_env %in% bfa_envs] <- df[[fit_col]][df$bfa_env %in% bfa_envs] / scale_factor
  return(df)
}


generate_resamples <- function(mu_vec, se_vec, bfa_envs, source, cluster, resamples = 100) {
  bfa_envs <- as.character(bfa_envs)
  var_matr <- matrix(rep(0, length(mu_vec)^2),
                     ncol = length(mu_vec))
  diag(var_matr) <- se_vec^2
  sim = paste0(c(1:resamples))
  df_sim <- data.frame(source,
                       cluster,
                       sim = sim,
                       MASS::mvrnorm(n = resamples,
                                     mu = mu_vec,
                                     Sigma = var_matr))
  names(df_sim) <- c("source", "cluster", "sim", bfa_envs)
  return(df_sim %>%
           tidyr::pivot_longer(cols = bfa_envs, names_to = "bfa_env", values_to = "s")
  )
}

plot_cluster_by_source <- function(df, source = "GlyEtOH", resamples = 100) {
  
  df %>%
    dplyr::filter(!is.na(cluster)) ->
    df
  
  df$source <- as.character(df$source)
  df$cluster <- as.character(df$cluster)
  df$bfa_env <- as.character(df$bfa_env)
  
  # take MVN samples from each vector of wmu_x and wse per source, cluster
  df %>%
    split(list(df$source, df$cluster), drop = TRUE) ->
    df_split
  
  df_split %>%
    lapply(function(x) generate_resamples(mu_vec = x[["wmu_x"]],
                                          se_vec = x[["wse"]],
                                          bfa_envs = x[["bfa_env"]],
                                          source = x[["source"]],
                                          cluster = x[["cluster"]],
                                          resamples = 50)) %>%
    do.call(rbind, .) ->
    df_resampled
  
  df %>%
    rescale_fitness(bfa_envs = c("FLC4", "CLM"), scale_factor = 4, fit_col = "wmu_x") %>%
    dplyr::filter(source == "GlyEtOH") ->
    df_rescaled
  
  df_resampled %>%
    dplyr::filter(source == "GlyEtOH") %>%
    rescale_fitness(bfa_envs = c("FLC4", "CLM"), scale_factor = 4, fit_col = "s") %>%
    ggplot(aes(x = bfa_env, y = s, group = sim)) +
    geom_hline(yintercept = 0, lty = "solid", col = "black", lwd = 0.25) +
    geom_vline(xintercept = "GlyEtOH", lty = "solid", col = "midnightblue", lwd = 5, alpha = 0.15) +
    #geom_hline(yintercept = c(-0.10, -0.05, 0.05, 0.10), lty = "solid", col = "gray30", lwd = 0.125) +
    geom_line(alpha = 0.2, col = "gray60") +
    geom_line(data = df_rescaled,
              aes(x = bfa_env, y = wmu_x, group = cluster), col = "darkorange2") +
    geom_linerange(data = df_rescaled,
                   aes(x = bfa_env, y = wmu_x, ymin = wmu_x - wse, ymax = wmu_x + wse, group = cluster),
                   col = "darkorange2") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA),
          #panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    facet_wrap(~ cluster, nrow = length(unique(df_resampled$cluster)))
  
}


# ==== #
# main #
# ==== #

"plot_by_cluster.R

Usage:
    plot_by_cluster.R [--help]
    plot_by_cluster.R [options] <fitness_file> <pleiotropy_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
Arguments:
    fitness_file                  fitness-by-cluster file (output 1 from summarise_clusters.R)
    pleiotropy_file               pleiotropy summary file (output 2 from summarise_clusters.R)
" -> doc

# define default args for debug_status == TRUE
args <- list(
  use_iva = TRUE,
  fitness_file = "data/combined/dBFA2_cutoff-5_compiled_data_by_barcode.csv",
  pleiotropy_file = "data/mutation_data/mutations_by_bc.csv",
  outdir = "output/figures/clusters"
)

debug_status <- FALSE

cat("\n*******************\n")
cat("* plot_by_cluster.R *\n")
cat("*********************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
