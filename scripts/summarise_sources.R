#! /usr/bin/env Rscript
# summarise_sources.R

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
    dplyr::select(-bfa_env, -s, -s_se) %>%
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
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid = element_blank()) +
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
  
  input_data <- read.table(arguments$input_file,
                           sep = ",",
                           header = T,
                           stringsAsFactors = F)
  
  input_data %>%
    summarise_sources() ->
    source_summaries
  
  bfa_prfx <- strsplit(basename(arguments$input_file), "_")[[1]][1]
  
  # fx these outputs
  source_summaries$source_summary %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$input_file),
              split_on = "_compiled",
              str_to_append = "_source_summaries_table")
  
  source_summaries$plot_table %>%
    write_out(out_dir = file.path(arguments$outdir, "tables"),
              base_name = basename(arguments$input_file),
              split_on = "_compiled",
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
  input_file = "data/combined/hBFA1_cutoff-5_compiled_data_by_barcode.csv",
  outdir = "output",
  gens = "8",
  exclude = "X48Hr"
)

debug_status <- FALSE

cat("\n***********************\n")
cat("* summarise_sources.R *\n")
cat("***********************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
