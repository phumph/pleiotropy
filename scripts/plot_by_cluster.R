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

theme_plt <- function(base_size = 10, base_family = "") {
  theme_minimal() %+replace%
    theme(axis.ticks = element_line(colour = "black", size = 0.2),
          panel.border = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "gray80", size = 0.25, linetype = "dashed"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8,
                                   family = "sans",
                                   color = "gray15"))
}

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggrepel)))

source(file.path("scripts/src/pleiotropy_functions.R"))

# ------------- #
# function defs #
# ------------- #

rescale_fitness <- function(df, bfa_envs, scale_factor, fit_col = "s", var_col, rescale_var = TRUE) {
  df[[fit_col]][df$bfa_env %in% bfa_envs] <- df[[fit_col]][df$bfa_env %in% bfa_envs] / scale_factor
  if (rescale_var == TRUE) {
    df[[var_col]][df$bfa_env %in% bfa_envs] <- df[[var_col]][df$bfa_env %in% bfa_envs] / scale_factor
  }
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


plot_means_by_source <- function(df) {
  
  source_order <- c("YPD", "SC",
                    "21C", "37C",
                    "02M_NaCl", "GlyEtOH",
                    "pH3_8", "pH7_3",
                    "FLC4", "CLM")
  
  df$source <- factor(df$source, levels = source_order)
  
  df %>%
    dplyr::filter(!source %in% c("FLC4", "CLM")) %>%
    ggplot(aes(x = wmu_x_home, y = wmu_x_away)) +
    facet_wrap(~ source, ncol = 2) +
    geom_hline(yintercept = 0, lwd = 0.5, linetype = "solid", color = "gray80") +
    geom_vline(xintercept = 0, lwd = 0.5, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(ymin = wmu_x_away - wse_away,
                      ymax = wmu_x_away + wse_away), col = "gray15", width = 0) +
    geom_errorbarh(aes(xmin = wmu_x_home - wse_home, 
                       xmax = wmu_x_home + wse_home), col = "gray15", height = 0) +
    geom_point(size = 2, col = "gray15") +
    coord_cartesian(xlim = c(0, 0.1),
                    ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = seq(-0.06, 0.06, 0.02)) +
    scale_x_continuous(breaks = seq(0, 0.1, 0.02)) +
    theme_plt() +
    #xlab("home fitness") +
    xlab("") +
    ylab("away fitness (mean)") +
    geom_label_repel(aes(label = cluster), label.size = NA, box.padding = 0,
                     fill = NA,
                     size = 3,
                     color = "gray40") ->
    mu_plot_no_drug
  
  df %>%
    dplyr::filter(source %in% c("FLC4", "CLM")) %>%
    ggplot(aes(x = wmu_x_home, y = wmu_x_away)) +
    facet_wrap(~ source, ncol = 2) +
    geom_hline(yintercept = 0, lwd = 0.5, linetype = "solid", color = "gray80") +
    geom_vline(xintercept = 0, lwd = 0.5, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(ymin = wmu_x_away - wse_away,
                      ymax = wmu_x_away + wse_away), col = "gray15", width = 0) +
    geom_errorbarh(aes(xmin = wmu_x_home - wse_home, 
                       xmax = wmu_x_home + wse_home), col = "gray15", height = 0) +
    geom_point(size = 2, col = "gray15") +
    coord_cartesian(xlim = c(0, 0.8),
                    ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = seq(-0.06, 0.06, 0.02)) +
    scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
    theme_plt() +
    xlab("home fitness") +
    ylab("") +
    geom_label_repel(aes(label = cluster), label.size = NA, box.padding = 0,
                     fill = NA,
                     size = 3,
                     color = "gray40") ->
    mu_plot_drug
  
  mu_plot <- ggpubr::ggarrange(plotlist = list(mu_plot_no_drug, mu_plot_drug),
                               nrow = 2, align = 'v',
                               heights = c(1, 0.31))
  return(mu_plot)
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


plot_pleio <- function(p_df, s_df, source = NULL) {
  
  #pal <- brewer.pal(name = "RdYlBu", n = 10)
  pal <- rev(pals::parula(n=11))
  
  # pleiotropy counts plot
  # ::::::::::::::::::::::::
  
  p_df %>%
    dplyr::select(source, cluster, n_pleio_pos, n_pleio_neg, n_pleio_net) %>%
    tidyr::pivot_longer(c(n_pleio_pos, n_pleio_neg, n_pleio_net),
                        names_to = "pleio_type", values_to = "pleio_score") %>%
    dplyr::mutate(source_cluster = paste0(source, "_", cluster)) ->
    p_df_pleio_long
  
  p_df_pleio_long %>%
    dplyr::filter(pleio_type == "n_pleio_net") %>%
    dplyr::arrange(source, pleio_score) %>%
    dplyr::select(source, source_cluster) ->
    clusters_for_order
  
  p_df_pleio_long$source_cluster <- factor(p_df_pleio_long$source_cluster,
                                           levels = clusters_for_order$source_cluster)
  
  hline_mids <- cumsum(table(clusters_for_order$source)) + 0.5
  
  fill_pal <- rev(pals::tol(n = length(unique(p_df_pleio_long$source))))
  #fill_pal <- rev(pals::kelly(n = length(unique(p_df_pleio_long$source))))
  
  p_df_pleio_long %>%
    dplyr::filter(pleio_type != "n_pleio_net") %>%
    ggplot(aes(x = pleio_type, y = source_cluster, size = pleio_score)) +
    geom_point(aes(fill = source), col = "black", pch = 21) +
    theme_plt() +
    theme(legend.position = "none",
          panel.border = element_blank()) +
    ylab("cluster ID") +
    xlab("") +
    scale_x_discrete(labels = c('-','+')) +
    scale_fill_manual(values = fill_pal) +
    geom_hline(yintercept = hline_mids, col = "black") ->
    pleio_plot_a
  
  # p_df_pleio_long %>%
  #   dplyr::filter(pleio_type == "n_pleio_net") %>%
  #   ggplot(aes(y = source_cluster, x = pleio_score, fill = source)) +
  #   #geom_bar(stat="identity", col = "gray15", lwd = 0.25) +
  #   geom_bar(stat="identity", col = "black", lwd = 0.25) +
  #   theme_plt() +
  #   theme(axis.text.y = element_blank(),
  #         legend.position = "none",
  #         panel.border = element_blank()) +
  #   # scale_fill_gradientn(colors = pal,
  #   #                      name="",
  #   #                      limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
  #   scale_fill_manual(values = fill_pal) +
  #   geom_vline(xintercept = 0) +
  #   ylab("") +
  #   xlab("net pleiotropy") +
  #   scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 2)) +
  #   geom_hline(yintercept = hline_mids, col = "black") ->
  #   pleio_plot_b

  p_df_pleio_long %>%
    dplyr::filter(pleio_type == "n_pleio_net") %>%
    ggplot(aes(y = source_cluster, x = pleio_score, fill = source)) +
    geom_vline(xintercept = 0) +
    #geom_bar(stat="identity", col = "gray15", lwd = 0.25) +
    #geom_bar(stat="identity", col = "black", lwd = 0.25) +
    geom_segment(aes(x = 0, xend = pleio_score, y = source_cluster, yend = source_cluster),
                 color = "gray15", lwd = 0.25) +
    geom_point(size = 2.5, pch = 21, col = "black") +
    theme_plt() +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank()) +
    # scale_fill_gradientn(colors = pal,
    #                      name="",
    #                      limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = fill_pal) +
    ylab("") +
    xlab("net pleiotropy") +
    scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 2)) +
    geom_hline(yintercept = hline_mids, col = "black") ->
    pleio_plot_b
  
  plots_b <- ggpubr::ggarrange(plotlist = list(pleio_plot_a, pleio_plot_b),
                               ncol = 2,
                               align = "h",
                               widths = c(0.6, 1))
  
  # fitness cluster plot
  # ::::::::::::::::::::::::
  
  # re-scale fitness for extreme environments
  s_df %>%
    rescale_fitness(bfa_envs = c("FLC4", "CLM"), scale_factor = 6,
                    fit_col = "wmu_x",
                    var_col = "wse") ->
    s_df_rescaled
  
  # join color for plotting clusters by pleiotropic net effect
  s_df_rescaled %>%
    dplyr::left_join(dplyr::select(p_df, source, cluster, n_pleio_net),
                     by = c("source", "cluster")) ->
    s_df_rescaled
  
  # re-order bfa-env
  bfa_levels <-  c("YPD", "SC",
                   "21C", "37C",
                   "02M_NaCl", "GlyEtOH",
                   "pH3_8", "pH7_3",
                   "FLC4", "CLM")
  s_df_rescaled$bfa_env <- factor(s_df_rescaled$bfa_env,
                                  levels = bfa_levels)
  
  s_df_rescaled %>%
    ggplot(aes(x = bfa_env, y = wmu_x, group = cluster, color = n_pleio_net)) +
    geom_hline(yintercept = 0, lty = "solid", col = "black", lwd = 0.25) +
    geom_vline(aes(xintercept = source), lty = "solid", col = "midnightblue", lwd = 5, alpha = 0.15) +
    geom_line() +
    geom_point(size = 1) +
    geom_linerange(aes(x = bfa_env, y = wmu_x, ymin = wmu_x - wse, ymax = wmu_x + wse)) +
    ylab("fitness change") +
    theme_plt() +
    facet_wrap(~ source, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_gradientn(colors = pal,
                          name="",
                          limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.04)) +
    coord_cartesian(ylim = c(-0.125, 0.125)) ->
    fit_plot
  
  return(list(
    fit_plot,
    plots_b
  ))
  # ggpubr::ggarrange(plotlist = list(fit_plot, plots_b), labels = c("a", "b"),
  #                   nrow = 2,
  #                   heights = pheights,
  #                   common.legend = TRUE, align = "h")
  # 
}


# ==== #
# main #
# ==== #


main <- function(arguments) {
  
  fit_df <- read.table(arguments$fitness_file,
                       sep = ",",
                       header = TRUE,
                       stringsAsFactors = FALSE)
  
  pleio_df <- read.table(arguments$pleiotropy_file,
                         sep = ",",
                         header = TRUE,
                         stringsAsFactors = FALSE)
  
  bfa_prfx <- strsplit(basename(arguments$fitness_file), "_")[[1]][1]
  
  strsplit(basename(arguments$fitness_file), "_")[[1]][1:2] %>% 
    paste0(collapse = "_") ->
    bfa_basename
  
  fig_out_base <- paste0(dirname(arguments$fitness_file),
                         "/",
                         bfa_basename,
                         "_") %>% gsub("tables", "figures", .)
  
  # mean fitness biplots
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  mu_plots_by_source <- plot_means_by_source(pleio_df)
  mu_plots_by_source %>%
    ggplot2::ggsave(filename = paste0(fig_out_base, "mu_plot.pdf"),
                    device = "pdf",
                    width = 4, height = 10, units = "in")
  
  # pleio plots by source
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  pleio_df %>%
    split(.$source) ->
    pleio_list
  
  fit_df %>%
    split(.$source) ->
    fit_list
  
  sources <- names(pleio_list)
  stopifnot(all(sources %in% names(fit_list)))
  pleio_plots <- vector("list", length = length(source_list))
  for (i in seq_along(sources)) {
    plot_filename <- file.path(arguments$outdir, 
                               paste0(bfa_prfx, "_", 
                                      sources[i],
                                      "_clust_fig.pdf"))
    pleio_plots[[i]] <- plot_pleio(p_df = pleio_list[[sources[i]]],
                                   s_df = fit_list[[sources[i]]],
                                   source = sources[i])
    
    pleio_plots[[i]][1] %>%
      ggplot2::ggsave(filename = plot_filename, device = "pdf",
                      width = 5
      )
    
  }
}



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
  fitness_file = "output/tables/dBFA2_cutoff-5_cluster_summaries_plot-data.csv",
  pleiotropy_file = "output/tables/dBFA2_cutoff-5_cluster_summaries_table.csv",
  outdir = "output/figures/clusters",
  png = TRUE
)

debug_status <- TRUE

cat("\n*******************\n")
cat("* plot_by_cluster.R *\n")
cat("*********************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
