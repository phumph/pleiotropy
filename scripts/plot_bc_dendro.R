#! /usr/bin/env Rscript

# plot_bc_dendro.R

# Inputs
# 1. compile data by barcode table for both ploidies
# 2. mutation data by barcode for both ploidies

# Outputs 
# 1. Heatmap sub-plots
# 2. Mutation heatmap sub-plot
# 3. Dendrogram sub-plots

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(ggdendro)))
suppressWarnings(suppressMessages(library(dendextend)))

source(file.path("scripts/src/pleiotropy_functions.R"))
source(file.path("scripts/src/plt_themes.R"))

# ------------- #
# function defs #
# ------------- #


full_rankify <- function(df) {
    df %>%
        tidyr::pivot_wider(id_cols = c(row_id, source, ploidy),
                           names_from = bfa_env,
                           values_from = s) %>%
        tidyr::pivot_longer(cols = names(.)[!names(.) %in% c("row_id", "source", "ploidy")],
                            names_to = "bfa_env",
                            values_to = "s") ->
        df_filled
    return(df_filled)
}


interpolate <- function(sub_df) {
    # Split by bfa_enf; mean-impute for NA values
    # Assumes each sub_df is of a single source
    stopifnot(length(unique(sub_df$source)) == 1)
    if (sum(is.na(sub_df$s)) == 0) {
        return(sub_df)
    }
    sub_sub_df <- split(sub_df, sub_df$bfa_env)
    lapply(sub_sub_df, function(x) {
        if (sum(is.na(x$s)) == 0) {
            return(x)
        }
        x$s[is.na(x$s)] <- mean(x$s, na.rm = TRUE)
        return(x)
    }) ->
        interpolated_envs
    return(
        do.call(rbind, interpolated_envs)
    )
}


interpolate_s <- function(df) {
    df_split <- split(df, df$source)
    df_split %>%
        lapply(interpolate) ->
        interpolated_srcs
    return(
        do.call(rbind, interpolated_srcs)
    )
}


apply_df_filters <- function(df, focal_ploidy = c("1N", "2N"),
                             remove_home = FALSE,
                             wgs_only = FALSE,
                             adapted_only = TRUE,
                             bfa_excludes = NULL) {
    if (remove_home == TRUE) {
        df %>%
            dplyr::mutate(is_home = ifelse(source == bfa_env, TRUE, FALSE)) %>%
            dplyr::filter(is_home == FALSE) ->
            df
    }
    if (wgs_only == TRUE) {
        df %>%
            dplyr::filter(has_wgs == TRUE) ->
            df
    }
    if (adapted_only == TRUE) {
        df %>%
            dplyr::filter(is_adapted == TRUE,
                          !is.na(cluster)) ->
            df
    }
    if (!is.null(bfa_excludes)) {
        df %>%
            dplyr::filter(!bfa_env %in% bfa_excludes) ->
            df
    }
    return(df %>% dplyr::filter(ploidy %in% focal_ploidy))
}


pre_cluster_for_dendro <- function(df, bfa_excludes = NULL) {
    # prep df into matrix for clustering
    if (!is.null(bfa_excludes)) {
        df %>%
            dplyr::filter(!bfa_env %in% bfa_excludes) ->
            df
    }
    df %>%
        dplyr::select(row_id, s, bfa_env) %>%
        tidyr::pivot_wider(id_cols = row_id,
                           names_from = bfa_env,
                           values_from = s) ->
        df_wide
    df_wide %>%
        dplyr::select(-row_id) %>%
        as.matrix ->
        mtr_wide
    row.names(mtr_wide) <- df_wide$row_id
    
    # do clustering
    df_dendro <- as.dendrogram(hclust(d = dist(x = mtr_wide)))
    env_dendro <- as.dendrogram(hclust(d = dist(x = t(mtr_wide))))
    env_dendro_dat <- dendro_data(env_dendro)
    dendro_dat <- dendro_data(df_dendro)
    # output factors in order for plotting
    return(
        list(
            "dends" = list(
                "bc_dendro" = df_dendro,
                "env_dendro" = env_dendro
            ),
            "names" = list(
                "bc_labels" = dendro_dat$labels$label,
                "bfa_labels" = env_dendro_dat$labels$label
            )
        )
    )
}


plot_bc_dendro <- function(dend, dend_labels, df, leaf_col_type) {
    make_dendro <- function(...) {
        dendro_plot <- plot(dend, horiz=TRUE)
        dendro_plot <- legend(x = 1.2, y = dim(source_labs)[1],
                              legend = unique(source_labs[, leaf_col_type]),
                              fill = unique(source_labs[, leaf_col]),
                              bty = "n",
                              cex = 0.5)
        #bar_colors <- ifelse(ploidies == "1N", "skyblue", "midnightblue")
        #env_colors <- source_labs$ploidy_col
        #colored_bars(colors = source_labs$source_col, dend = dend, rowLabels = c("ploidy"))
    }
    
    source_cols <- pals::tol(n = length(unique(df$source)))
    source_col_df <- data.frame(source = unique(df$source),
                                source_col = source_cols,
                                stringsAsFactors = FALSE)
    
    data.frame(row_id = dend_labels) %>%
        dplyr::left_join(dplyr::select(df, row_id, source, ploidy) %>%
                             unique()) %>%
        dplyr::left_join(source_col_df) %>%
        dplyr::mutate(ploidy_col = ifelse(ploidy == "1N", "skyblue", "midnightblue")) ->
        source_labs
    leaf_col = ifelse(leaf_col_type == "source", "source_col", "ploidy_col")
    leaf_colors = source_labs[, leaf_col]
    suppressWarnings(
        dend %>%
            set("leaves_pch", 19) %>% 
            set("leaves_cex", 1) %>% 
            set("leaves_col", leaf_colors) %>%
            set("labels", "") ->
            dend
    )
    return(make_dendro)
}


rescale_s <- function(df) {
    df %>%
        tidyr::pivot_wider(id_cols = c(row_id, source, ploidy),
                           names_from = bfa_env,
                           values_from = s) ->
        df_wide
    
    df_wide[ , -c(1,2,3)] %>%
        as.matrix() %>%
        apply(2, scale) %>%
        data.frame() ->
        df_wide_scaled
    
    df_wide_scaled$row_id <- df_wide$row_id
    df_wide_scaled %>%
        tidyr::pivot_longer(cols = names(df_wide_scaled)[names(df_wide_scaled) != "row_id"],
                            names_to = "bfa_env",
                            values_to = "s") %>%
        dplyr::mutate(bfa_env = gsub("^X", "", bfa_env)) ->
        df_long
    return(df_long)
}


plot_heatmap <- function(df) {
    df %>%
        ggplot(aes(x = bfa_env, y = row_id, fill = s)) +
        geom_tile() +
        theme_plt() +
        theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_blank(),
              legend.position = "top",
              axis.ticks = element_blank()) +
        scale_fill_gradientn(colors = pals::brewer.puor(256),
                             limits = c(-3,3),
                             na.value = "gray60") +
        ylab("") +
        xlab("") ->
        heatmap_plot
    return(heatmap_plot)
}


prep_mutation_data <- function(df, bc_labels) {
    df %>%
        dplyr::group_by(GENE) %>%
        dplyr::summarise(nhits = length(GENE)) %>%
        dplyr::filter(nhits > 2) %>%
        dplyr::arrange(nhits) ->
        multi_hit_genes
    
    missing_rows <- data.frame(
        row_id = paste0(
            bc_labels[
                !bc_labels %in% df$row_id
            ]
        )
    )
    df %>%
        dplyr::mutate(row_id = Full.BC) %>%
        dplyr::group_by(row_id, GENE) %>%
        dplyr::summarise(hits = length(GENE)) %>%
        dplyr::filter(GENE %in% multi_hit_genes$GENE) %>%
        tidyr::pivot_wider(id_cols = row_id, names_from = GENE, values_from = hits, values_fill = 0) ->
        mdf_tmp
    
    mdf_tmp %>%
        dplyr::bind_rows(missing_rows) %>%
        tidyr::pivot_longer(cols = multi_hit_genes$GENE, names_to = "GENE", values_to = "hits") %>%
        dplyr::filter(row_id %in% bc_labels) ->
        mut_hits
    
    mut_hits$row_id <- factor(mut_hits$row_id, levels = bc_labels)
    mut_hits$GENE <- factor(mut_hits$GENE, levels = rev(paste0(multi_hit_genes$GENE)))
    mut_hits %>%
        dplyr::mutate(plot_hit = ifelse(hits > 0, hits, NA)) ->
        mut_hits
    return(mut_hits)
}


make_mutation_plot <- function(df) {
    df %>%
        dplyr::mutate(plot_hit_binary = ifelse(!is.na(plot_hit), 0.1, NA)) %>%
        ggplot(aes(x = GENE, y = row_id)) +
        geom_hline(yintercept = df$row_id[!is.na(df$plot_hit)], lwd = 0.25, alpha = 0.5) +
        geom_point(aes(size = plot_hit), alpha = 0.5) +
        theme_plt() +
        theme(axis.text.y = element_blank(),
              legend.position = "none",
              panel.grid.major.y = element_blank(),
              axis.ticks.y = element_blank()) +
        theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5)) +
        ylab("") +
        xlab("") ->
        mut_heatmap
    return(mut_heatmap)
}


main <- function(arguments) {
    
    # Prepare fitness data
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # load fitness files
    fit_dat <- load_files_from_arg(arguments$fitness_files)
    
    # add ploidy flag
    arguments$fitness_files %>%
        stringr::str_extract_all(stringr::regex(".BFA")) %>%
        unlist() ->
        assays_fit
    ploidies = c("hBFA" = "1N", "dBFA" = "2N")
    ploidy_dat <- ploidies[assays_fit]
    purrr::map2(fit_dat, ploidy_dat, function(x, y) dplyr::mutate(x, ploidy = y)) ->
        fit_dat
    
    fit_dat %>%
        dplyr::bind_rows() ->
        fit_df_full
    
    # transform fitness data in preparation for dendro
    fit_df_full %>%
        apply_df_filters(focal_ploidy = ploidies,
                         remove_home = FALSE,
                         wgs_only = FALSE,
                         adapted_only = TRUE,
                         bfa_excludes = "37C_Stan") ->
        fit_df_filt
    
    fit_df_filt %>%
        dplyr::mutate(row_id = Full.BC) %>%
        dplyr::select(row_id, source, bfa_env, s, ploidy) %>%
        dplyr::arrange(row_id) ->
        df_tmp
    
    # interpolate NA fitness values as mean of source
    # within a given bfa_env
    df_tmp %>%
        full_rankify() %>%
        interpolate_s() ->
        fit_df_interpolated
    
    # perform clustering on barcodes and output list with
    # 1. dendro object
    # 2. ordered factor levels (BC) for heatmap
    dendro_and_labels <- pre_cluster_for_dendro(df = fit_df_interpolated)
    
    dendro_and_labels$dends$bc_dendro %>%
        plot_bc_dendro(dend_labels = dendro_and_labels$names$bc_labels,
                       df = df_tmp,
                       leaf_col_type = "source") -> 
        plot_dendro_source
    dendro_and_labels$dends$bc_dendro %>%
        plot_bc_dendro(dend_labels = dendro_and_labels$names$bc_labels,
                       df = df_tmp,
                       leaf_col_type = "ploidy") -> 
        plot_dendro_ploidy
    
    # save dendro sub-plots
    cowplot::plot_grid(plot_dendro_source) %>% 
        ggsave(filename = file.path(args$outdir,"dendro_by_bc_source.pdf"),
               width = 6,
               height = 12,
               units = "in",
               device = "pdf"
        )
    cowplot::plot_grid(plot_dendro_ploidy) %>% 
        ggsave(filename = file.path(args$outdir, "dendro_by_bc_ploidy.pdf"),
               width = 6,
               height = 12,
               units = "in",
               device = "pdf"
        )
    
    # Produce fitness heatmap
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::

    # re-scale fitness per column (mean centered, sigma-scaled)
    df_tmp %>%
        rescale_s() ->
        df_scaled

    df_scaled$row_id <- factor(df_scaled$row_id, levels = dendro_and_labels$names$bc_labels)
    df_scaled$bfa_env <- factor(df_scaled$bfa_env, levels = dendro_and_labels$names$bfa_labels)
    heatmap_plot <- plot_heatmap(df_scaled)
    heatmap_plot %>%
        ggsave(filename = file.path(args$outdir, "heatmap_by_bc.pdf"),
               width = 1.65,
               height = 12,
               units = "in",
               device = "pdf"
        )
    
    # Prepare mutation data
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    mut_dat <- read.csv(arguments$mutation_file)
    mut_dat %>%
        prep_mutation_data(bc_labels = dendro_and_labels$names$bc_labels) ->
        mut_df_for_plot
    
    # Produce mutation heatmap
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::

    mut_df_for_plot %>%
        make_mutation_plot() ->
        mut_plot
    
    ggpubr::ggarrange(plotlist = list(heatmap_plot, mut_plot),
                      widths = c(1, 8), ncol = 2,
                      align = "hv", common.legend = TRUE) %>%
        ggsave(filename = file.path(args$outdir, "heatmap_with_muts.pdf"),
               height = 12,
               width = 18,
               units = "in",
               device = "pdf")
}


# ==== #
# main #
# ==== #


"plot_bc_dendro.R

Usage:
    plot_bc_dendro.R [--help]
    plot_bc_dendro.R [options] <fitness_files> <mutation_file>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
Arguments:
    fitness_files                 space-separated quoted string of fitness files
    mutation_file                 mutations_by_bc file
" -> doc

# define default args for debug_status == TRUE
args <- list(
    fitness_files = "data/combined/hBFA1_cutoff-5_compiled_data_by_barcode.csv data/combined/dBFA2_cutoff-5_compiled_data_by_barcode.csv",
    mutation_file = "data/mutation_data/mutations_by_bc.csv",
    outdir = "output/figures",
    png = FALSE
)

debug_status <- FALSE

cat("\n********************\n")
cat("* plot_bc_dendro.R *\n")
cat("********************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
