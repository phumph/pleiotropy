# pleiotropy_functions.R
#
# Imports for pleiotropy analysis.

# ======== #
# fun defs #
# ======== #

run_args_parse <- function(arguments = NULL, debug_status = FALSE) {
  if (debug_status == TRUE) {
    stopifnot(!is.null(arguments) & is.list(arguments))
  } else if (debug_status == FALSE) {
    arguments <- docopt(doc)
  }
  return(arguments)
}


OpenRead <- function(arg) {
  if (arg %in% c("-", "/dev/stdin")) {
    file("stdin", open = "r")
  } else {
    file(arg, open = "r")
  }
}


prep_fit_matrix <- function(dat,
                            is_neutral = FALSE,
                            is_autodip = FALSE,
                            means_only = FALSE,
                            neutral_col,
                            autodip_tag,
                            excludes,
                            iva_s,
                            gens) {

  gens <- as.double(gens)
  
  # take input counts file and generate list of 2 matrixes:
    # (n x k) matrix of k-variate mean vectors for n barcodes
    # (n x k) k-variate sigma vectors for n barcodes
  avg_type <- ifelse(iva_s == TRUE, 'iva_s', 'ave')
  fit_col  <- ifelse(iva_s == TRUE, 'iva_s', 'ave_s')
  err_col  <- ifelse(iva_s == TRUE, 'iva_s_err', 'ave_err')

  focal_cols <- c(grep('\\.BC',   names(dat), value = T),
                  grep('Subpool', names(dat), value = T),
                  grep('\\.R[0-9]', names(dat), invert = T, value = T),
                  grep(paste0(avg_type),  names(dat), value = T)) %>% unique()

  # focal_cols <- focal_cols[!grepl(paste0(unlist(excludes), collapse='|'), focal_cols)]
  focal_cols <- focal_cols[!grepl(excludes, focal_cols)]
  dat <- dat[ , names(dat) %in% focal_cols]

  # examine list of neutrals if flag present
  if (is_neutral == TRUE && !is.null(neutral_col)) {
    # filter to neutrals
    dat <- dat[dat$Subpool.Environment == neutral_col, ]
  }

  # examine list of neutrals if flag present
  if (is_autodip == TRUE && !is.null(autodip_tag)) {
    # filter to autodips
    dat <- dat[grepl(autodip_tag, dat$Which.Subpools), ]
  }

  row.names(dat) <- dat$Full.BC

  # grab relevant fitness columns from input
  df_mat <-
    dat %>%
    dplyr::select(grep(paste0(fit_col,'$|',err_col,'$'), names(dat), value = T)) %>%
    as.matrix()

  df_mat[is.infinite(df_mat)] <- NA
  df_mat <- df_mat[complete.cases(df_mat), ]

  # split into means and sigmas:
  if (means_only == FALSE) {
    df_list <- list(means  = df_mat[, grep(paste0(fit_col,'$'), colnames(df_mat))] / gens,
                    sigmas = df_mat[, grep(paste0(err_col,'$'), colnames(df_mat))] / gens)
    return(df_list)
  } else if (means_only == TRUE) {
    means <- df_mat[, grep(paste0(fit_col,'$'), colnames(df_mat))] / gens
    return(means)
  }
}


filter_neutral_outliers <- function(dfs,
                                    reps_iter,
                                    cutoff,
                                    init_cutoff = 0.05,
                                    GENS = 8,
                                    R = 500,
                                    outlier_p_cutoff = 0.05) {

  # for each barcode left out:
  # generate predictive distribution by sampling from vector of means and stderrs

  means  <- dfs[["means"]]
  sigmas <- dfs[["sigmas"]]

  # find distance of each barcode to mean; sort by desc(distance)
  starting_dists <- mahalanobis(means,
                                colMeans(means),
                                cov(means))
  means <-
    means %>%
    cbind(starting_dists)

  sigmas <-
    sigmas %>%
    cbind(starting_dists)

  means  <- means[order(means[, 'starting_dists'], decreasing = TRUE), ]
  sigmas <- sigmas[order(sigmas[, 'starting_dists'], decreasing = TRUE), ]

  # determine quick and dirty the set I'm going to evaluate
  # all with Chisq p-value of < init_cutoff
  init_pvals <- pchisq(means[, 'starting_dists'], df = dim(means)[2] - 1, lower.tail = FALSE)

  # grab barcodes with p-values <= init_cutoff:
  candidate_outliers <- names(init_pvals[init_pvals <= init_cutoff])

  cat(sprintf("%s initial outlier candidates detected.\n", length(candidate_outliers)))
  cat("Progressively eliminating outlier barcodes...\n")

  pvals <- NULL
  removed_outliers <- NULL

  # go through each candidate outlier
  # remove it
  # generate more robust testing distribution
  # flag to remove given p-value
  for (i in seq_along(candidate_outliers)) {

    # remove focal barcode
    means2 <- means[!row.names(means) == candidate_outliers[i],
                    grep('starting_dists', colnames(means), invert = TRUE)]

    sigmas2 <- sigmas[!row.names(sigmas) == candidate_outliers[i],
                      grep('starting_dists', colnames(sigmas), invert = TRUE)]

    test_bc <- means[row.names(means) == candidate_outliers[i],
                     grep('starting_dists', colnames(means), invert = TRUE)]

    # generate multivariate normal random variables for each vector:
    res_mat <- NULL

    # sample MVN random vectors from remaining neutrals
    for (ii in seq_along(row.names(means2))) {
      res_mat <- rbind(res_mat, MASS::mvrnorm(n = R, means2[ii, ], diag(sigmas2[ii, ]^2)))
    }

    # calculate distance of each barcode in set
    D <- mahalanobis(rbind(test_bc,res_mat),
                     colMeans(res_mat),
                     cov(res_mat))

    # generate distance of test barcode
    D_test <- mahalanobis(test_bc,
                          colMeans(res_mat),
                          cov(res_mat))

    # calculate p-value of test_bc dist
    P_test <- sum(D >= D_test) / length(D)

    # grow vector of P-values for records
    pvals <- append(pvals, P_test)

    # if P_test below threshold, remove from list
    if (P_test <= outlier_p_cutoff) {

      cat(sprintf('\tEliminated %s\n', candidate_outliers[i]))

      removed_outliers <- append(removed_outliers, candidate_outliers[i])

      means <- means[!row.names(means) == candidate_outliers[i],
                     grep('starting_dists', colnames(means), invert = TRUE)]

      sigmas <- sigmas[!row.names(sigmas) == candidate_outliers[i],
                        grep('starting_dists', colnames(sigmas), invert = TRUE)]
    }
  }

  cat(sprintf('Eliminated %s total barcodes flagged as outliers.\n', length(removed_outliers)))

  final_neuts <- row.names(means)[!row.names(means) %in% removed_outliers]

  return(final_neuts)
}


generate_neutral_test_distn <- function(dfs, neutrals, reps_final) {

  means  <- dfs[["means"]]
  sigmas <- dfs[["sigmas"]]

  # subset to neutral list
  means2  <-  means[row.names(means) %in% neutrals, ]
  sigmas2 <- sigmas[row.names(sigmas) %in% neutrals, ]

  # generative model of distribution of k-variate vectors
  res_mat <- NULL
  for (ii in seq_along(row.names(means2))) {
    res_mat <- rbind(res_mat, MASS::mvrnorm(n = reps_final, means2[ii, ], diag(sigmas2[ii, ]^2)))
  }

  # generate testing distribution of distances to apply to all barcodes
  D <- mahalanobis(res_mat,
                   colMeans(res_mat),
                   cov(res_mat))

  L1 <- list(col_means = colMeans(res_mat),
             cov_matr  = cov(res_mat),
             distances = D)

  return(L1)
}


flag_adapteds <- function(matr, testing_distn, cutoff) {

  # takes mean vector of all barcodes, neutral test distn,
  # and compares distance of each barcode to the neutral set.
  # for now, does not handle missing values but simply excludes them.
  # I'll need to update the script to do this soon.

  # calculate distance of each barcode versus neutral set:
  D <- mahalanobis(matr,
                   testing_distn$col_means,
                   testing_distn$cov_matr)

  # for each: determine p-value within testing distribution
  p_vals <- sapply(D, function(x) sum(testing_distn$distances >= x) / length(testing_distn$distances))
  is_adapted <- D > quantile(testing_distn$distances, probs = 1 - cutoff)

  return(data.frame(Full.BC = row.names(matr),
                    matr,
                    maha_dist = D,
                    dist_pval = p_vals,
                    is_adapted))
}


filter_to_focal_bcs <- function(input_df,
                                adapt_col   = 'is_adapted',
                                neutral_col = 'neutral_set',
                                autodip_col = 'autodip',
                                retain_neutrals = FALSE,
                                retain_autodips = TRUE,
                                retain_adapteds = TRUE) {
  
  stopifnot(all(c(adapt_col, neutral_col) %in% names(input_df)))
  stopifnot(sum(retain_neutrals, retain_autodips, retain_adapteds) > 0)
  
  # filter out autodips if present
  if ((retain_autodips == FALSE) & (autodip_col %in% names(input_df))) {
    input_df %>%
      dplyr::filter(!!sym(autodip_col) == FALSE) ->
      input_df
  }
  # assemble set based on neutral and adapted criteria
  if (retain_adapteds == TRUE & retain_neutrals == TRUE) {
    input_df %>%
      dplyr::filter(!!sym(adapt_col) == TRUE) %>%
      dplyr::filter(!!sym(neutral_col) == TRUE) ->
      input_df
    return(input_df)
  }
  if (retain_adapteds == TRUE & retain_neutrals == FALSE) {
    input_df %>%
      dplyr::filter(!!sym(adapt_col) == TRUE) %>%
      dplyr::filter(!!sym(neutral_col) == FALSE) ->
      input_df
    return(input_df)
  }
  if (retain_adapteds == FALSE & retain_neutrals == TRUE) {
    input_df %>%
      dplyr::filter(!!sym(adapt_col) == FALSE) %>%
      dplyr::filter(!!sym(neutral_col) == TRUE) ->
      input_df
    return(input_df)
  }
  return(input_df)
}
