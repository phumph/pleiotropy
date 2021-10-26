#! /usr/bin/env Rscript

# plot_cluster_dendro.R

# Inputs
# 1. cluster summaries result table for both ploidies
# 2. compiled mutation data for both ploidies

# Outputs 
# 1. individual plots per cluster per source environment
# 2. summary plot of pleiotropic effects per source environment

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggrepel)))
library(ggdendro)
library(dendextend)

source(file.path("scripts/src/pleiotropy_functions.R"))
source(file.path("scripts/src/plt_themes.R"))

# ------------- #
# function defs #
# ------------- #

plot_dendro <- function(df, mut_df, focal_ploidy = "2N", remove_home = FALSE) {
  
  if (remove_home == TRUE) {
    df %>%
      dplyr::filter(is_home == FALSE) ->
      df
  }
  
  df %>%
    dplyr::filter(ploidy %in% focal_ploidy) %>%
    dplyr::mutate(row_id = paste0(source, "_", ploidy, "_", cluster)) %>%
    dplyr::select(row_id, source, ploidy, bfa_env, wmu_x) %>%
    dplyr::arrange(row_id) ->
    df_labs
  
  df_labs %>%
    dplyr::select(-source, -ploidy) %>%
    tidyr::pivot_wider(id_cols = row_id, names_from = bfa_env, values_from = wmu_x, ) ->
    df_wide
  
  # for clustering, remove columns with NA:
  na_cols <- colSums(apply(df_wide, 2, is.na))
  df_wide %>%
    dplyr::select(-row_id, -names(na_cols[na_cols > 0])) ->
    df_wide_pruned
  
  df_matr <- as.matrix(df_wide_pruned)
  row.names(df_matr) <- df_wide$row_id
  
  df_dendro <- as.dendrogram(hclust(d = dist(x = df_matr)))
  env_dendro <- as.dendrogram(hclust(d = dist(x = t(df_matr))))
  env_dendro_dat <- dendro_data(env_dendro)
  dendro_dat <- dendro_data(df_dendro)
  dendro_dat$labels$label <- as.character(dendro_dat$labels$label)
  
  source_cols <- pals::tol(n = length(unique(df_labs$source)))
  source_col_df <- data.frame(source = unique(df_labs$source),
                              col = source_cols,
                              stringsAsFactors = FALSE)

  dendro_dat$labels %>% 
    dplyr::select(label) %>%
    dplyr::left_join(dplyr::select(df_labs, row_id, source) %>% 
                       unique(), by = c("label" = "row_id")) %>%
    dplyr::left_join(source_col_df) ->
    source_labs

  df_dendro %>%
    set("leaves_pch", 19) %>% 
    set("leaves_col", paste0(source_labs$col)) %>%
    set("labels", "") ->
    #set("branches_col", source_labs$color)
    dend
  
  if (all(c("1N", "2N") %in% focal_ploidy)) {
    ploidy_labels <- ifelse(grepl("_1N_", dendro_dat$labels$label), "black", "skyblue")
    colored_bars(colors = ploidy_labels,
                 dend = dend,
                 horiz = TRUE,
                 y_shift = 0.025,
                 rowLabels = "")  
  }

  df_wide[ , -1] %>%
    apply(2, scale) %>%
    data.frame() ->
    df_wide_scaled

  df_wide_scaled$row_id <- df_wide$row_id

  # plot heatmap
  df_wide_scaled %>%
    dplyr::select(-names(na_cols[na_cols > 0])) %>%
    tidyr::pivot_longer(cols = names(df_wide_scaled)[names(df_wide_scaled) != "row_id"],
                        names_to = "bfa_env",
                        values_to = "s") %>%
    dplyr::mutate(bfa_env = gsub("^X", "", bfa_env)) ->
    df_long
  
  df_long$row_id <- factor(df_long$row_id, levels = dendro_dat$labels$label)
  df_long$bfa_env <- factor(df_long$bfa_env, levels = paste0(env_dendro_dat$labels$label))
  df_long %>%
    ggplot(aes(x = bfa_env, y = row_id, fill = s)) +
    geom_tile() +
    theme_plt() +
    theme(axis.text.x = element_text(sze = 7, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank(),
          legend.position = "top") +
    scale_fill_gradientn(colors = pals::brewer.puor(256)) +
    ylab("") +
    xlab("") ->
    heatmap_plot
  
  # MUTATION PLOT
  # :::::::::::::::::::::::::::::::::::::::::::::::::
  
  # prepare mutation data
  mdf %>%
    dplyr::filter(PLOIDY %in% focal_ploidy) %>%
    dplyr::group_by(GENE) %>%
    dplyr::summarise(nhits = length(GENE)) %>%
    dplyr::filter(nhits > 2) %>%
    dplyr::arrange(nhits) ->
    multi_hit_genes
  
  # add missing row_ids
  missing_rows <- data.frame(row_id = paste0(dendro_dat$labels$label[!dendro_dat$labels$label %in% mdf_2$row_id]))
  
  mdf %>%
    dplyr::filter(PLOIDY %in% focal_ploidy) %>%
    dplyr::mutate(row_id = paste0(source, "_", PLOIDY, "_", cluster)) %>%
    dplyr::group_by(row_id, GENE) %>%
    dplyr::summarise(hits = length(GENE)) %>%
    dplyr::filter(GENE %in% multi_hit_genes$GENE) %>%
    tidyr::pivot_wider(id_cols = row_id, names_from = GENE, values_from = hits, values_fill = 0) ->
    mdf_tmp
  
  mdf_tmp %>%
    dplyr::bind_rows(missing_rows) %>%
    tidyr::pivot_longer(cols = multi_hit_genes$GENE, names_to = "GENE", values_to = "hits") ->
    mut_hits
  
  # re-order row_id by dendro
  mut_hits$row_id <- factor(mut_hits$row_id, levels = levels(df_long$row_id))
  
  # re-order GENE by total count
  mut_hits$GENE <- factor(mut_hits$GENE, levels = rev(paste0(multi_hit_genes$GENE)))
  
  mut_hits %>%
    dplyr::mutate(plot_hit = ifelse(hits > 0, hits, NA)) ->
    mut_hits
  
  mut_hits %>%
    ggplot(aes(x = GENE, y = row_id)) +
    geom_hline(yintercept = mut_hits$row_id[!is.na(mut_hits$plot_hit)], lwd = 0.25, alpha = 0.5) +
    geom_point(aes(size = plot_hit), alpha = 0.5) +
    theme_plt() +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_blank()) +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5)) +
    #scale_fill_gradientn(colors = pals::brewer.reds(256)) +
    ylab("") +
    xlab("") ->
    mut_heatmap
  
  # MAKE DENDROGRAM PLOT
  # :::::::::::::::::::::::::::::::::::::::::::::::::
  
  make_dendro <- function() {
    dendro_plot <- plot(dend, horiz=TRUE)
    dendro_plot <- legend(x = 0.9, y = dim(source_labs)[1]/1.25, #"topleft",
                          legend = unique(source_labs$source),
                          fill = unique(source_labs$col),
                          bty = "n",
                          cex = 0.5)
  }
  
  ggpubr::ggarrange(plotlist = list(heatmap_plot, mut_heatmap),
                    widths = c(1, 3.5),
                    align = "hv", common.legend = TRUE) %>%
    ggsave(filename = "output/figures/dBFA2_heatmaps.pdf",
           height = plot_params$dBFA_map_w,
           width = plot_params$dBFA_map_h,
           units = "in",
           device = "pdf")
  
  # save dendro
  cowplot::plot_grid(make_dendro) %>%
    ggsave(filename = "output/figures/dBFA2_dendro.pdf",
           width = plot_params$dBFA_dendro_w,
           height = plot_params$dBFA_dendro_h,
           units = "in",
           device = "pdf")
}

plot_params <- list(
  dBFA_map_w = 5,
  dBFA_map_h = 6,
  dBFA_dendro_w = 5,
  dBFA_dendro_h = 7,
  hBFA_map_w = 5,
  hBFA_map_h = 6,
  hBFA_dendro_w = 5,
  hBFA_dendro_h = 7,
)


main <- function(arguments) {
  
  # Prepare input data
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  # load fitness files
  fit_dat <- load_files_from_arg(arguments$fitness_files)
  
  # load mutation files
  mut_dat <- load_files_from_arg(arguments$mutation_files)
  
  # add ploidy flag
  arguments$fitness_files %>%
    stringr::str_extract_all(stringr::regex(".BFA")) %>%
    unlist() ->
    assays_fit
  
  arguments$mutation_files %>%
    stringr::str_extract_all(stringr::regex(".BFA")) %>%
    unlist() ->
    assays_mut

  stopifnot(all(assays_mut == assays_fit))

  ploidies = c("hBFA" = "1N", "dBFA" = "2N")
  ploidy_dat <- ploidies[assays_fit]

  purrr::map2(fit_dat, ploidy_dat, function(x, y) dplyr::mutate(x, ploidy = y)) ->
    fit_dat

  fit_dat %>%
    dplyr::bind_rows() ->
    fit_df

  mut_dat %>%
    dplyr::bind_rows() ->
    mut_df

  # Cluster clusters
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::

  plot_dendro(df = fit_df, mdf = mut_df, focal_ploidy = "2N", remove_home = TRUE)
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
    -p --png                      Flag to determine output device as PNG [default: PDF]
Arguments:
    fitness_file                  fitness-by-cluster file (output 1 from summarise_clusters.R)
    pleiotropy_file               pleiotropy summary file (output 2 from summarise_clusters.R)
" -> doc

# define default args for debug_status == TRUE
args <- list(
  fitness_files = "output/tables/hBFA1_cutoff-5_cluster_summaries_plot-data.csv output/tables/dBFA2_cutoff-5_cluster_summaries_plot-data.csv",
  mutation_files = "data/mutation_data/hBFA1_mutations_by_cluster.csv data/mutation_data/dBFA2_mutations_by_cluster.csv",
  outdir = "output/figures",
  png = FALSE
)

debug_status <- TRUE

cat("\n*************************\n")
cat("* plot_cluster_dendro.R *\n")
cat("*************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
