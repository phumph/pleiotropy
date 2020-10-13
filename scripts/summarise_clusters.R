#! /usr/bin/env Rscript

# summarise_cluster.R

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

source(file.path("scripts/src/pleiotropy_functions.R"))


# ------------- #
# function defs #
# ------------- #

weighted_mean <- function(x, se, stderr = FALSE) {
  if (length(x) == 1) {
    return(c("wmu_x" = x, "wse" = se))  
  }
  w = 1 / se^2
  wmu_x = sum(x*w) / sum(w)
  n = length(x)
  n_eff = sum(w)^2 / sum(w^2)
  wvar = (sum(w * (x - wmu_x)^2) / sum(w)) * (n_eff / (n_eff - 1))
  wse = sqrt(wvar)
  if (stderr == TRUE) {
    return(c("wmu_x" = wmu_x, "wse" = wse / sqrt(n_eff)))
  }
  return(c("wmu_x" = wmu_x, "wse" = wse))
}


summarise_clusters <- function(df) {

  # calculate weighted mean and var of each cluster
  df %>%
    dplyr::filter(is_adapted == TRUE, !is.na(cluster)) ->
    df_filt  
  
  # calculate num. BCs per cluster
  df_filt %>%
    dplyr::group_by(source, cluster) %>%
    dplyr::summarise(n_bcs = length(unique(Full.BC))) ->
    bcs_per_cluster
  bcs_per_cluster$cluster <- as.character(bcs_per_cluster$cluster)
  
  # apply calculation to each level of source, cluster, bfa_env
  df_split <- split(df_filt, list(df_filt$source, df_filt$cluster, df_filt$bfa_env))
  lapply(df_split, function(x) weighted.mean(x = x[["s"]], se = x[["s_se"]])) %>%
    do.call(rbind, .) -> df_filt_2
  
  # reconstitute results into data.frame
  df_filt_2 <- df_filt_2[complete.cases(df_filt_2),]
  strsplit(row.names(df_filt_2), "\\.") %>% 
    do.call(rbind, .) %>% 
    as.data.frame() ->
    df_cols
  names(df_cols) <- c("source", "cluster", "bfa_env")
  df_res <- data.frame(df_cols, df_filt_2, stringsAsFactors = FALSE)
  
  # count up nominally pleiotropic effects
  z_cutoff <- 1.96
  df_res$z = df_res$wmu_x / df_res$wse
  df_res$pleio <- 0
  df_res$pleio[df_res$z > z_cutoff] <- 1
  df_res$pleio[df_res$z < -z_cutoff] <- -1
  
  # calculate average fitness accross away environments
  df_res$is_home <- paste0(df_res$source) == paste0(df_res$bfa_env)
  
  df_res %>%
    dplyr::filter(is_home == FALSE) -> df_tmp
  
  df_tmp_split <- split(df_tmp, list(df_tmp$source, df_tmp$cluster))
  lapply(df_tmp_split, function(x) weighted_mean(x = x[["wmu_x"]], se = x[["wse"]], stderr = TRUE)) %>%
    do.call(rbind, .) -> df_tmp_res
  
  df_tmp_res <- df_tmp_res[complete.cases(df_tmp_res),]
  strsplit(row.names(df_tmp_res), "\\.") %>% 
    do.call(rbind, .) %>% 
    as.data.frame() ->
    df_cols
  names(df_cols) <- c("source", "cluster")
  df_tmp_res <- data.frame(df_cols, df_tmp_res, stringsAsFactors = FALSE)
  df_tmp_res %>%
    dplyr::rename(wmu_x_away = wmu_x, wse_away = wse) ->
    df_fit_away_res
  
  # summarise cluster info per source
  df_res %>%
    dplyr::group_by(source, cluster) %>%
    dplyr::summarise(n_pleio_pos = sum(pleio == 1),
                     n_pleio_neg = sum(pleio == -1),
                     n_pleio_none = sum(pleio == 0),
                     n_pleio_any = sum(abs(pleio)),
                     n_pleio_net = sum(pleio),
                     n_bfa_envs = n()) ->
    df_pleio
  
  df_res %>%
    dplyr::filter(is_home == TRUE) %>%
    dplyr::select(source, cluster, wmu_x, wse) %>%
    dplyr::rename(wmu_x_home = wmu_x,
                  wse_home = wse) ->
    df_fit_home_res
  
  suppressWarnings(
    df_pleio %>%
      dplyr::left_join(df_fit_home_res, by = c("source", "cluster")) %>%
      dplyr::left_join(df_fit_away_res, by = c("source", "cluster")) %>%
      dplyr::left_join(bcs_per_cluster, by = c("source", "cluster")) ->
      df_res_byclust
  )
  
  stopifnot(all(df_pleio$n_pleio_pos + df_pleio$n_pleio_neg + df_pleio$n_pleio_none == df_pleio$n_bfa_envs))
  
  return(df_res_byclust)
}


main <- function(arguments) {

  input_data <- read.table(arguments$input_file,
                           sep = ",",
                           header = T,
                           stringsAsFactors = F)
  
  input_data %>%
    summarise_clusters() ->
    cluster_summaries

  bfa_prfx <- strsplit(basename(arguments$input_file), "_")[[1]][1]

  cluster_summaries %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$input_file),
              str_to_append = "_cluster_summaries_table")

  if (grepl("dBFA2", bfa_prfx)) {
    plot_width <- 7
    plot_height <- 5
  } else {
    plot_width <- 5
    plot_height <- 4
  }

  cluster_summaries$plot %>%
    ggsave(filename = file.path(arguments$outdir,
                                "figures",
                                paste0(bfa_prfx,
                                       "_cluster_summaries_plot",
                                       ".pdf")),
           width = plot_width,
           height = plot_height,
           device = "pdf",
           units = "in")

}

# ==== #
# main #
# ==== #

"summarise_clusters.R

Usage:
    summarise_clusters.R [--help]
    summarise_clusters.R [options] <input_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations

Arguments:
    input_file                    Input data file (combined fitness, cluster, and mutations data)
" -> doc

# define default args for debug_status == TRUE
args <- list(
  use_iva = TRUE,
  input_file = "data/combined/dBFA2_cutoff-5_compiled_data_by_barcode.csv",
  outdir = "output",
  gens = "8",
  exclude = "X48Hr"
)

debug_status <- TRUE

cat("\n**********************\n")
cat("* summarise_csources.R *\n")
cat("************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
