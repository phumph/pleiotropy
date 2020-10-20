#! /usr/bin/env Rscript

# summarise_cluster.R

# This script takes input fitness data (post-clustering)
# as well as mutation data
# and generates cluster-level statistics 
# about fitness variation.
#
# Outputs 
# 1. bfa-level data for spaghetti plots
# 2. cluster-level data summarised for plotting
# 3. cluster-and-mutation links for assessing mutual information


# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))

source(file.path("scripts/src/pleiotropy_functions.R"))

Z_CUTOFF<- 1.96
BFA_EXCL <- "37C_Stan"

# ------------- #
# function defs #
# ------------- #

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
  df_split <- split(df_filt, list(df_filt$source, df_filt$cluster, df_filt$bfa_env), drop = TRUE)
  lapply(df_split, function(x) weighted_mean(x = x[["s"]], se = x[["s_se"]])) %>%
    do.call(rbind, .) -> df_filt_2
  
  # reconstitute results into data.frame
  strsplit(row.names(df_filt_2), "\\.") %>% 
    do.call(rbind, .) %>% 
    as.data.frame() ->
    df_cols
  names(df_cols) <- c("source", "cluster", "bfa_env")
  df_res <- data.frame(df_cols, df_filt_2, stringsAsFactors = FALSE)
  
  # count up nominally pleiotropic effects
  z_cutoff <- Z_CUTOFF
  df_res$z = df_res$wmu_x / df_res$wse
  df_res$pleio <- 0
  df_res$pleio[df_res$z > z_cutoff] <- 1
  df_res$pleio[df_res$z < -z_cutoff] <- -1
  
  # calculate average fitness accross away environments
  df_res$is_home <- paste0(df_res$source) == paste0(df_res$bfa_env)
  
  df_res %>%
    dplyr::filter(is_home == FALSE) -> 
    df_tmp
  
  df_tmp_split <- split(df_tmp, list(df_tmp$source, df_tmp$cluster), drop = TRUE)
  lapply(df_tmp_split, function(x) weighted_mean(x = x[["wmu_x"]], se = x[["wse"]], stderr = FALSE)) %>%
    do.call(rbind, .) -> df_tmp_res
  
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
  
  stopifnot(all(df_pleio$n_pleio_pos + df_pleio$n_pleio_neg + df_pleio$n_pleio_none == df_pleio$n_bfa_envs))
  
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
  
  df_res_byclust %>%
    dplyr::mutate(wsd_away = wse_away,
                  wse_away = wse_away / sqrt(n_bfa_envs)) ->
    df_res_byclust
  
  return(list(cluster_summary = df_res_byclust,
              bfa_summary = df_res))
}


parse_muts_for_output <- function(mut_df, fit_df, assay) {
  
  if (grepl("dBFA", assay)) {
    ploidy <- "2N"
  } else if (grepl("hBFA", assay)) {
    ploidy <- "1N"
  } else {
    ploidy <- NA
  }
  
  source_meta <- fit_df[, c("Full.BC", "source", "cluster")] %>% unique()
  focal_cols <- c("Full.BC", "PLOIDY", "VARIANT_NAME",
                  "AMINO_ACID_CHANGE", "REF", "ALT",
                  "GENE", "INDEL", "FUNCTIONAL_CATEGORY",
                  "EFFECT", "BIOTYPE", "IMPACT")
  
  mut_df %>%
    dplyr::rename(VARIANT_NAME = Tag,
                  FUNCTIONAL_CATEGORY = FUNCLASS) %>%
    dplyr::mutate(PLOIDY = ploidy,
                  INDEL = ifelse(nchar(REF) != nchar(ALT), TRUE, FALSE)) %>%
    dplyr::arrange(Full.BC) ->
    mut_parsed_df
  
  mut_parsed_df <- mut_parsed_df[, names(mut_parsed_df) %in% focal_cols]
    
  # join with source and cluster data
  mut_parsed_df %>%
    dplyr::left_join(source_meta, by = c("Full.BC")) ->
    mut_w_meta_df
  
  mut_w_meta_df %>%
    dplyr::filter(!is.na(source), !is.na(cluster)) ->
    mut_w_meta_filt_df
  
  return(mut_w_meta_filt_df)
}


main <- function(arguments) {
  
  input_data <- read.table(arguments$input_file,
                           sep = ",",
                           header = TRUE,
                           stringsAsFactors = FALSE)
  
  
  
  bfa_prfx <- strsplit(basename(arguments$input_file), "_")[[1]][1]
  
  mutation_data %>%
    parse_muts_for_output(fit_df = input_data, assay = bfa_prfx) ->
    mutations_parsed
  
  input_data %>%
    dplyr::filter(!grepl(paste0(arguments$exclude, collapse = "|"), source),
                  !grepl(paste0(BFA_EXCL, collapse = "|"), bfa_env)) ->
    input_data
  
  input_data %>%
    summarise_clusters() ->
    cluster_summaries
  
  cluster_summaries[[2]] %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$input_file),
              str_to_append = "_cluster_summaries_table")
}

# ==== #
# main #
# ==== #

"summarise_clusters.R

Usage:
    summarise_clusters.R [--help]
    summarise_clusters.R [options] <combined_file> <mutations_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations
Arguments:
    input_file                    Input fitness data file (combined fitness and cluster data)
    mutations_file                Mutations data file (mutations listed by barcode)
" -> doc

# define default args for debug_status == TRUE
args <- list(
  use_iva = TRUE,
  input_file = "data/combined/dBFA2_cutoff-5_compiled_data_by_barcode.csv",
  mutations_file = "data/mutation_data/mutations_by_bc.csv",
  outdir = "output",
  exclude = "48Hr"
)

debug_status <- TRUE

cat("\n**********************\n")
cat("* summarise_clusters.R *\n")
cat("************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
