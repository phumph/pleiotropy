#!/usr/bin/env Rscript
# calc_pw_clust_dist.R

# ------ #
# header #
# ------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(furrr)))

source(file.path("scripts/src/pleiotropy_functions.R"))

# will fail if run within Rstudio
# comment out the line below if interactively running code:
future::plan(multiprocess)


# ------------- #
# function defs #
# ------------- #

x <- rnorm(8, 0, 1)
y <- rnorm(8, 0, 1)
d <- dist_euclidean(x, y)


dist_euclidean <- function(x, y) {
    return(
        sqrt(sum((x - y)^2))
    )
}


prep_df_for_calc <- function(df){
    # gather wmu and wse fields as wide-format Dfs
    df %>%
        dplyr::select(-wse) %>%
        tidyr::pivot_wider(names_from = bfa_env, values_from = wmu_x) ->
        df_wmu
    df %>%
        dplyr::select(-wmu_x) %>%
        tidyr::pivot_wider(names_from = bfa_env, values_from = wse) ->
        df_wse
    
    # remove any bfa_envs with missing data
    cols_with_missing <- names(which(colSums(is.na(df_wmu[, -1])) > 0))
    df_wmu %>%
        dplyr::select(-cols_with_missing) ->
        df_wmu
    df_wse %>%
        dplyr::select(-cols_with_missing) ->
        df_wse
    
    # generate lists in preparation for map functions
    l_wmu <- split(df_wmu[, -1], df_wmu$clust_id)
    l_wse <- split(df_wse[, -1], df_wmu$clust_id)
    l_id <- split(df_wmu$clust_id, df_wmu$clust_id)
    
    # generate list of pairs of wmu and wse list for each cluster
    purrr::map2(l_wmu, l_wse, function(x, y) {
        return(list("wmu" = as.vector(t(x)), "wse" = as.vector(t(y))))
    }) -> L
    return(L)
}


make_pairwise_list <- function(list_) {
    # generate cartesian product of each single list element
    LL <- purrr::cross2(list_, list_)
    # same as above but with list element names
    names_ <- purrr::cross2(names(list_), names(list_))
    # generate comparison string to name cartesian product list
    LL_names <- purrr::map(names_, function(x) {
        return(paste0(x[[1]], ":", x[[2]]))
    })
    names(LL) <- LL_names
    return(LL)
}


calc_avg_diff <- function(l_, dist_fun = dist_euclidean, N = 50) {
    stopifnot(length(l_) == 2)
    # generate replicate mvn random vectors for each pairwise combination of list elements
    purrr::map(l_, function(x) {
        replicate(N, purrr::map2(x$wmu, x$wse, function(wmu, wse) { rnorm(1, wmu, wse) }) %>%
                      unlist())
    }) ->
        rmvn_vectors
    # calculate pairwise distances for each pair of vectors
    dists <- apply(rmvn_vectors[[1]], 2, dist_fun, rmvn_vectors[[2]])
    dist_mu <- mean(dists)
    dist_sdev <- sd(dists)
    return(c("mu" = dist_mu, "sigma" = dist_sdev))
}


parse_names <- function(df) {
    row.names(df) %>%
        strsplit(":") %>%
        purrr::map(sort) %>%
        do.call(rbind, .) %>%
        as.data.frame() ->
        names_df
    names(names_df) <- c("x_comp", "y_comp")
    names_df %>%
        dplyr::select(x_comp) %>%
        tidyr::separate(x_comp, c("x_ploidy", "x_source_cluster"),
                        sep = "_",
                        remove = FALSE,
                        extra = "merge") %>%
        tidyr::separate(x_source_cluster, "x_source", sep = "_", remove = FALSE, extra = "drop") %>%
        dplyr::mutate(x_cluster = strsplit(x_source_cluster, "_") %>% unlist() %>% last) %>%
        dplyr::select(-x_source_cluster) ->
        x_names
    names_df %>%
        dplyr::select(y_comp) %>%
        tidyr::separate(y_comp, c("y_ploidy", "y_source_cluster"),
                        sep = "_",
                        remove = FALSE,
                        extra = "merge") %>%
        tidyr::separate(y_source_cluster, "y_source", sep = "_", remove = FALSE, extra = "drop") %>%
        dplyr::mutate(y_cluster = strsplit(y_source_cluster, "_") %>% unlist() %>% last) %>%
        dplyr::select(-y_source_cluster) ->
        y_names
    df_parsed <- dplyr::bind_cols(x_names, y_names) %>% unique()
    df_parsed %>%
        dplyr::mutate(join_key = paste0(x_comp, ":", y_comp)) ->
        df_parsed
    df %>%
        dplyr::mutate(join_key = row.names(.)) ->
        df
    df_parsed %>%
        dplyr::left_join(df, by = "join_key") %>%
        dplyr::select(-join_key) ->
        df_final
    return(df_final)
}


annotate_comparisons <- function(df) {
    df %>%
        dplyr::mutate(is_self = x_comp == y_comp,
                      same_env = x_source == y_source,
                      ploidy = ifelse(x_ploidy != y_ploidy, "1N_2N", paste0(x_ploidy, "_", y_ploidy))) ->
        df
    df$ploidy <- factor(df$ploidy, levels = c("1N_1N", "2N_2N", "1N_2N"))
    return(df)
}


plot_distances <- function(df) {
    sources = c("GlyEtOH", "FLC4", "CLM", "21C")
    df %>%
        dplyr::filter(#is_self == FALSE,
                      same_env == TRUE,
                      x_source %in% sources) ->
        df_filt
    
    df_filt %>%
        ggplot() +
        #geom_histogram(aes(x = mu, fill = ploidy), bins = 50, alpha = 0.8) +
        geom_density(aes(x = mu, fill = ploidy), alpha = 0.9, col = NA) +
        facet_grid(ploidy ~ x_source) +
        theme_bw() +
        theme(strip.background = element_blank(),
              legend.position = "none",
              panel.grid.minor = element_blank()) +
        scale_fill_manual(values = pals::brewer.set2(3)) +
        xlab("Euclidean dist between clusters") ->
        dist_plot
    return(dist_plot)
}


main <- function(arguments) {
    # process input files
    fit_dat <- load_files_from_arg(arguments$fitness_files)
    
    # gather ploidy info from filenames
    arguments$fitness_files %>%
        stringr::str_extract_all(stringr::regex(".BFA")) %>%
        unlist() ->
        assays_fit
    ploidies = c("hBFA" = "1N", "dBFA" = "2N")
    ploidy_dat <- ploidies[assays_fit]
    purrr::map2(fit_dat, ploidy_dat, function(x, y) dplyr::mutate(x, ploidy = y)) ->
        fit_dat
    
    # gather two ploidy Dfs into one
    fit_dat %>%
        dplyr::bind_rows() ->
        fit_df_full
    
    # reduce Df shape and generate composite unique row ID
    fit_df_full %>%
        dplyr::mutate(clust_id = paste0(ploidy, "_", source, "_", cluster)) %>%
        dplyr::select(-source, -cluster, -ploidy, -pleio, -is_home, -z) ->
        df_slim
    
    # generate paired wmu and wse list for each cluster_id
    list_by_cluster <- prep_df_for_calc(df_slim)
    
    # generate list containing all pairwise combinations of cluster lists
    pairwise_list <- make_pairwise_list(list_by_cluster)
    
    # generate pairwise distances using list as input
    distances_list <- purrr::map(pairwise_list, calc_avg_diff, dist_euclidean)
    names(distances_list) <- names(pairwise_list)
    do.call(rbind, distances_list) %>% 
        as.data.frame() -> 
        res_df
    # reconcile names of comparitors; eliminate lower triangle
    # and retain only upper and diag entries
    res_df %>%
        parse_names() %>%
        annotate_comparisons() ->
        res_df_parsed
    
    res_df_parsed %>%
        plot_distances() ->
        dist_plot
    
    dist_plot %>%
        ggsave(filename = file.path(arguments$outdir, "pairwise_cluster_dist_plot.png"),
                                    device = "png", 
                                    width = 6,
                                    height = 4,
                                    dpi = 300)
}


# ==== #
# main #
# ==== #

"calc_pw_clust_dist.R

Usage:
    calc_pw_clust_dist.R [--help]
    calc_pw_clust_dist.R [options] <fitness_files>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
Arguments:
    fitness_files                 space-separated quoted string of combined fitness files
" -> doc

# define default args for debug_status == TRUE
args <- list(
    fitness_files = "output/tables/hBFA1_cutoff-5_cluster_summaries_plot-data.csv output/tables/dBFA2_cutoff-5_cluster_summaries_plot-data.csv",
    outdir = "output/figures",
    png = TRUE
)

debug_status <- TRUE

cat("\n************************\n")
cat("* calc_pw_clust_dist.R *\n")
cat("************************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")

