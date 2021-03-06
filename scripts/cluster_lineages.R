#! /usr/bin/env Rscript

# cluster_lineages.R

# This script takes input fitnesses and
# generates cluters for each source-ploidy lineage set.
#
# The script also calculates the cluster-wise weighted (e.g., Maha) mean and se.
# and outputs this as a second file

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(docopt)))

source(file.path("scripts/src/pleiotropy_functions.R"))

# ------------------------- #
# function definitions      #
# ------------------------- #

get_neutral_tree_depth <- function(df, method = "euclidean", quant = 0.95) {

  df %>%
    dist(method = method) %>%
    hclust() ->
    neut_clust

  q <- stats::quantile(neut_clust$height, prob = quant)

  return(as.numeric(q))
}


get_clusters_by_cut <- function(df = NULL, method = "euclidean", cut_at = NULL) {

  stopifnot(!is.null(cut_at) & !is.null(df))

  df <- df[complete.cases(df), ]

  df %>%
    dist(method = method) %>%
    hclust() %>%
    cutree(h = cut_at) ->
    cuts

  return(data.frame(Full.BC = row.names(df),
                    cluster = cuts,
                    stringsAsFactors = FALSE))
}


get_neutral_cov_matr <- function(mat) {
  cov(mat)
}


maha_mean <- function(x, x_var, covm) {

  if (all(is.vector(x), is.vector(x_var))) {
    sigma_hat <- covm
    diag(sigma_hat) <- x_var

    return(list(mean = x,
                covm = sigma_hat))
  }

  for (i in seq_along(nrow(x))) {
    sigma_i <- covm
    diag(sigma_i) <- x_var[i, ]
    sigma_i_inv <- solve(sigma_i)
    maha_x <- sigma_i_inv %*% x[i, ]

    if (i == 1) {
      sigma_hat_inv <- sigma_i_inv
      maha_x_tot <- maha_x
    } else {
      sigma_hat_inv <- sigma_hat_inv + sigma_i_inv
      maha_x_tot <- maha_x_tot + maha_x
    }
    sigma_hat <- solve(sigma_hat_inv)
    maha_mean_vec <- sigma_hat %*% maha_x_tot
    maha_mean_vec2 <- c(maha_mean_vec)
    names(maha_mean_vec2) <- row.names(maha_mean_vec)
    maha_mean_vec <- maha_mean_vec2
  }

  return(list(mean = maha_mean_vec,
              covm = sigma_hat))
}


calc_clust_dist <- function(...) {

  args <- list(...)

  stopifnot(all(c("C_i", "C_j", "dist_fun", "covm") %in% names(args)))
  stopifnot(exists(args$dist_fun))

  dist_fun <- get(args$dist_fun)

  x_i <- dist_fun(x = args$C_i$x_i, x_var = args$C_i$x_i_var, covm = args$covm)
  x_j <- dist_fun(x = args$C_j$x_j, x_var = args$C_j$x_j_var, covm = args$covm)

  d_ij = t(x_i$mean - x_j$mean) %*%
    solve(x_i$covm + x_j$covm) %*%
    (x_i$mean - x_j$mean)

  # equivalent to mahalanobis(x_i$mean, x_j$mean, cov = x_i$covm + x_j$covm)

  return(d_ij)
}


numerize <- function(x) {
  if (nrow(x) == 1){
    return(as.vector(x))
  } else {
    return(as.matrix(x))
  }
}

reduce_clusters_by_dist <- function(clust_vec = NULL,
                                    mean_matr = NULL,
                                    var_matr = NULL,
                                    covm = NULL,
                                    se = FALSE,
                                    dist_fun = "maha_mean",
                                    p_cutoff = 0.99) {

  stopifnot(all(!is.null(covm), !is.null(mean_matr), !is.null(var_matr)))

  df <- ncol(mean_matr)
  dist_cutoff <- qchisq(p = p_cutoff, df = df)

  if (is.null(clust_vec)) {
    clust_vec <- seq_len(dim(mean_matr)[1])
  } else {
    stopifnot(is.vector(clust_vec) &
                length(clust_vec) == nrow(mean_matr) &
                nrow(mean_matr) == nrow(var_matr))
  }

  if (se == TRUE) {
    var_matr <- var_matr^2
  }

  while (length(unique(clust_vec)) > 1) {

    DISTS <- list()

    clust_ids <- unique(clust_vec)

    clust_vec %>%
      plyr::mapvalues(from = clust_ids,
                      to = seq_len(length(unique(clust_vec)))) ->
      clust_vec

    n_clusts <- length(unique(clust_vec))
    clust_ids <- unique(clust_vec)

    cat(sprintf("\tNum. clusters = %s\n", n_clusts))

    # calculate pairwise dist between each cluster
    for (i in 1:(n_clusts - 1)) {
      if (i %in% clust_vec == FALSE) {
        i <- i + 1
      }
      c_i <- clust_ids[i]
      x_i <- as.matrix(mean_matr[clust_vec == c_i, ])
      x_i_var <- as.matrix(var_matr[clust_vec == c_i, ])
      C_i <- list(x_i = x_i,
                  x_i_var = x_i_var)
      for (j in (i + 1):n_clusts) {
        if (j %in% clust_vec == FALSE) {
          j <- j + 1
        }
        c_j <- clust_ids[j]
        x_j <- as.matrix(mean_matr[clust_vec == c_j, ])
        x_j_var <- as.matrix(var_matr[clust_vec == c_j, ])
        C_j <- list(x_j = x_j,
                    x_j_var = x_j_var)
        # calc dist
        d_ij = calc_clust_dist(C_i = C_i,
                               C_j = C_j,
                               covm = covm,
                               dist_fun = "maha_mean")
        DISTS <- append(DISTS, list(data.frame(i = c_i, j = c_j, d_ij = d_ij)))
      }
    }
    D_res <- do.call(rbind, DISTS)
    D_min <- min(D_res$d_ij)
    if (D_min >= dist_cutoff) {
      break
    }
    row_min <- which(D_res$d_ij == D_min)
    clust_vec[clust_vec %in% c(D_res$i[row_min],
                               D_res$j[row_min])] <- D_res$i[row_min]
  }
  cat("Done!\n")
  if (all(!is.null(row.names(mean_matr)), !is.null(row.names(var_matr)))) {
    res <- data.frame(Full.BC = row.names(mean_matr),
                      cluster = clust_vec,
                      stringsAsFactors = FALSE)
    return(res)
  }
  return(clust_vec)
}


average_clusters <- function(clusters, mean_matr, var_matr, covm) {

  if (is.data.frame(clusters) & "cluster" %in% names(clusters)) {
    clust_ids <- clusters$cluster
  } else {
    clust_ids <- clusters
  }

  df_full <- data.frame()

  for (clust in seq_along(unique(clust_ids))) {

    clust_mean <- maha_mean(x = as.matrix(mean_matr[clust_ids == clust, ]),
                            x_var = as.matrix(var_matr[clust_ids == clust, ]),
                            covm = covm)

    names(clust_mean$mean) %>%
      sapply(function(x) {
        x %>%
          strsplit("\\.") %>%
          unlist() %>%
          dplyr::first()
      }) -> envs

    names(clust_mean$mean) <- envs

    data.frame(bfa_env = names(clust_mean$mean),
               s = clust_mean$mean,
               s_se = sqrt(diag(clust_mean$covm)),
               cluster = clust,
               n_bcs = length(clust_ids[clust_ids == clust])) ->
      df_tmp

    df_full %>%
      dplyr::bind_rows(df_tmp) ->
      df_full
  }
  return(df_full)
}


# ------------------------- #
# main def                  #
# ------------------------- #

main <- function(arguments) {

  # check + grab input files
  if (sum(!sapply(arguments$infiles, file.exists)) > 0) {
    stop(
        "One or more infile does not exist. Please check infile path and try again.",
         call. = FALSE
       )
  }

  # read infiles
  infile <- read.table(arguments$infile,
                       sep = ",",
                       header = T,
                       stringsAsFactors = F)

  # grab barcodes to retain after filtering:
  adapted_df <-
    infile %>%
    dplyr::filter(Subpool.Environment != "not_read") %>%
    filter_to_focal_bcs(
      retain_neutrals = FALSE,
      retain_adapteds = TRUE,
      retain_autodips = FALSE
    )

  neutral_df <-
    infile %>%
    filter_to_focal_bcs(
      retain_neutrals = TRUE,
      retain_adapteds = FALSE,
      retain_autodips = FALSE
    )

  # find tree depth of neutrals with mahalanobis distance:
  neutral_df %>%
    prep_fit_matrix(means_only = FALSE,
                    excludes = arguments$exclude,
                    iva_s = arguments$use_iva,
                    gens = arguments$gens) ->
    neutral_df_matr

  neutral_df_matr$means %>%
    get_neutral_tree_depth(quant = 0.95) ->
    neutral_clust_height

  neutral_df_matr$means %>%
    get_neutral_cov_matr() ->
    neutral_cov

  # use neutral tree depth
  # to define initial clusters for each environment
  adapted_df %>%
    split(adapted_df$Subpool.Environment) ->
    adapted_by_env

  adapted_df_w_clust <- list()
  clust_wise_means <- list()

  for (env in seq_along(adapted_by_env)) {

    cat(sprintf("\nWorking on env %s...\n", names(adapted_by_env)[env]))

    adapted_by_env[[env]] %>%
      prep_fit_matrix(means_only = FALSE,
                      excludes = arguments$exclude,
                      iva_s = arguments$use_iva,
                      gens = arguments$gens) ->
      matr_prepped
    
    if (all(lapply(matr_prepped, nrow) == 1)) {
      clusts_final <- data.frame(Full.BC = row.names(matr_prepped$means),
                                 cluster = "1",
                                 stringsAsFactors = FALSE)
    } else {
      matr_prepped$means %>%
        get_clusters_by_cut(cut_at = neutral_clust_height) ->
        clusts_initial
      
      clusts_initial$cluster %>%
        reduce_clusters_by_dist(mean_matr = matr_prepped$means,
                                var_matr = matr_prepped$sigmas^2,
                                covm = neutral_cov) ->
        clusts_final
    }

    adapted_by_env[[env]] %>%
      dplyr::left_join(clusts_final, by = "Full.BC") ->
      adapted_df_w_clust[[env]]
    
    clusts_final %>%
      average_clusters(mean_matr = matr_prepped$means,
                       var_matr = matr_prepped$sigmas^2,
                       covm = neutral_cov) ->
      clust_wise_mean
    
    clust_wise_mean$source <- names(adapted_by_env)[env]
    clust_wise_means[[env]] <- clust_wise_mean
  }

  adapted_df_w_clust_full <- do.call(rbind, adapted_df_w_clust)
  clust_wise_means_full   <- do.call(rbind, clust_wise_means)

  adapted_df_w_clust_full %>%
    write_out(out_dir = arguments$outdir,
              base_name = basename(arguments$infile),
              str_to_append = "_adapted_w_clusts")

  clust_wise_means_full %>%
    write_out(out_dir = arguments$outdir,
              base_name = basename(arguments$infile),
              str_to_append = "_adapted_w_clust_means")
}


# ------------------------- #
# main                      #
# ------------------------- #

"cluster_lineages.R

Usage:
    cluster_lineages.R [--help]
    cluster_lineages.R [options] <infile>

Options:
    -h --help                     Show this screen.
    -o --outdir=<outdir>          Output directory [default: ./]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations

Arguments:
    infile                        Input file(s) containing fitness calls for BFA run.
" -> doc

# define default args for debug_status == TRUE
arguments <- list(
  use_iva = TRUE,
  infile = "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_autodips.csv",
  outdir = "data/fitness_data/fitness_calls",
  gens = "8",
  exclude = "CLM|FLC4|Stan|48Hr"
)

debug_status <- FALSE
arguments <- run_args_parse(arguments, debug_status)

cat("\n**********************\n")
cat("* cluster_lineages.R *\n")
cat("**********************\n\n")

main(arguments)

cat("**Script completed successfully!**\n\n")
