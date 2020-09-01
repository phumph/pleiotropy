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

source(file.path("scripts/src/pleiotropy_functions.R"))

# ------------- #
# function defs #
# ------------- #


summarise_sources <- function(df) {

  df %>%
    dplyr::select(-bfa_env, -s) %>%
    unique() %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n_bcs = length(unique(Full.BC)),
                     n_bcs_adapted = sum(is_adapted),
                     p_bcs_adapted = round(n_bcs_adapted / n_bcs, 2),
                     n_wgs = sum(has_wgs),
                     p_wgs = round(n_wgs / n_bcs, 2),
                     n_bcs_adapted_wgs = sum(is_adapted == TRUE & has_wgs == TRUE),
                     p_adapted_w_wgs = round(n_bcs_adapted_wgs / n_bcs_adapted, 2),
                     n_clusters = length(unique(cluster)))  %>%
    dplyr::arrange(desc(n_bcs_adapted)) ->
    df2

  # fit summary plot
  df2 %>%
    dplyr::select(-n_wgs, -p_bcs_adapted) %>%
    tidyr::pivot_longer(cols = c("n_bcs", "n_bcs_adapted", "n_bcs_adapted_wgs"),
                        names_to = "bc_set",
                        values_to = "bc_count") ->
    df3

  # re-order factors:
  df3 %>%
    dplyr::filter(bc_set == "n_bcs_adapted") %>%
    dplyr::arrange(desc(bc_count)) %>%
    dplyr::select(source) ->
    source_order

  df3$source <- factor(df3$source, levels = (source_order$source))

  sources_plot <-
    ggplot() +
    geom_bar(data = df3, aes(x = source, y = bc_count, fill = bc_set),
             stat = "identity",
             position = "dodge", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = c("gray40", "dodgerblue", "darkorange2"),
                      name = "bc_set") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("num. BCs") +
    xlab("source env") +
    annotate(geom = "text",
             label = paste0(df2$p_bcs_adapted),
             x = paste0(df2$source), y = -12, size = 2) +
    annotate(geom = "text",
             label = paste0(df2$p_adapted_w_wgs),
             x = paste0(df2$source), y = -26, size = 2)

  return(list(source_summary = df2,
              plot = sources_plot,
              plot_table = df3))
}


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
    prep_fit_matrix(means_only = TRUE,
                    excludes = arguments$exclude,
                    iva_s = arguments$use_iva,
                    gens = arguments$gens) %>%
    as.data.frame() %>%
    dplyr::mutate(Full.BC = row.names(.)) %>%
    tidyr::pivot_longer(cols = grep("_s$", names(.), value = T),
                        names_to = "bfa_env", values_to = "s") %>%
    dplyr::left_join(
      fitness_data[, c("Full.BC", "Subpool.Environment", fit_dat_cols)],
      by = "Full.BC") ->
    fitness_matr_long

  fitness_matr_long %>%
    normalize_envs() ->
    fitness_matr_norm

  # exclude autodiploitds if hBFA
  if (grepl("hBFA", basename(arguments$fitness_file))) {
    fitness_matr_norm %>%
      dplyr::filter(autodip == FALSE) %>%
      dplyr::select(-autodip) ->
      fitness_matr_norm
  }

  fitness_matr_norm %>%
    summarise_sources() ->
    source_summaries

  bfa_prfx <- strsplit(basename(arguments$fitness_file), "_adapteds")[[1]][1]

  source_summaries$source_summary %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$fitness_file),
              str_to_append = "_source_summaries_table")

  source_summaries$plot_table %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$fitness_file),
              str_to_append = "_source_summaries_plot-data")

  if (grepl("dBFA2", bfa_prfx)) {
    plot_width <- 7
    plot_height <- 5
  } else {
    plot_width <- 5
    plot_height <- 4
  }

  source_summaries$plot %>%
    ggsave(filename = file.path(arguments$outdir,
                                "figures",
                                paste0(bfa_prfx,
                                       "_source_summaries_plot",
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
    summarise_clusters.R [options] <fitness_file> <cluster_file> <mutations_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
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
  fitness_file = "data/fitness_data/fitness_calls/dBFA2_cutoff-5_adapteds.csv",
  cluster_file = "data/fitness_data/fitness_calls/dBFA2_cutoff-5_adapted_w_clusts.csv",
  mutations_file = "data/mutation_data/mutations_by_bc.csv",
  outdir = "output",
  gens = "8",
  exclude = "X48Hr"
)

debug_status <- FALSE

cat("\n**********************\n")
cat("* summarise_csources.R *\n")
cat("************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
