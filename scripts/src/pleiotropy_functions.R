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


prep_fit_matrix <- function(dat,
                            is_neutral = FALSE,
                            is_autodip = FALSE,
                            means_only = FALSE,
                            neutral_col,
                            autodip_tag,
                            excludes = NULL,
                            iva_s,
                            gens) {
  
  gens <- as.double(gens)
  
  # take input counts file and generate list of 2 matrixes:
  # (n x k) k-variate mean vectors for n barcodes
  # (n x k) k-variate sigma vectors for n barcodes
  
  avg_type <- ifelse(iva_s == TRUE, 'iva_s', 'ave')
  fit_col  <- ifelse(iva_s == TRUE, 'iva_s', 'ave_s')
  err_col  <- ifelse(iva_s == TRUE, 'iva_s_err', 'ave_err')
  
  focal_cols <- c(grep('\\.BC',   names(dat), value = T),
                  grep('Subpool', names(dat), value = T),
                  grep('\\.R[0-9]', names(dat), invert = T, value = T),
                  grep(paste0(avg_type),  names(dat), value = T)) %>% unique()
  
  if (!is.null(excludes)) {
    focal_cols <- focal_cols[!grepl(excludes, focal_cols)]
  }
  
  dat <- dat[, names(dat) %in% focal_cols]
  
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
    dplyr::select(grep(paste0(fit_col,'$|',err_col,'$'), names(dat), value = T))
  
  df_mat <- sapply(df_mat, function(x) replace(x, is.infinite(x), NA), simplify = FALSE) %>%
    data.frame()
  
  row.names(df_mat) <- row.names(dat)
  
  if (is.data.frame(df_mat) && dim(df_mat)[1] > 1) {
    df_mat <- df_mat[complete.cases(df_mat), ]
  } else {
    stopifnot(sum(is.na(df_mat)) == 0)
  }
  
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
                                    cutoff = 0.05,
                                    init_cutoff = 0.1,
                                    GENS = 8,
                                    R = 1000) {
  
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
    if (P_test <= cutoff) {
      
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
  
  # calculate distance of each barcode versus neutral set:
  D <- mahalanobis(matr,
                   testing_distn$col_means,
                   testing_distn$cov_matr)
  
  # for each: determine p-value within testing distribution
  test_distn_len <- length(testing_distn$distances)
  p_vals <- sapply(D, function(x) {
    sum(testing_distn$distances > x) / test_distn_len
  })
  
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


normalize_envs <- function(df,
                           source_col = "Subpool.Environment",
                           remove_str_from_source = "_alpha$",
                           remove_str_from_bfa_env = "^X") {
  
  stopifnot(is.data.frame(df) & all(c(source_col, "bfa_env") %in% names(df)))
  
  s_str <- stringr::str_extract(df$bfa_env, "\\..+") %>%
    unique()
  
  stopifnot(length(s_str) == 1)
  
  if (!is.na(s_str)) {
    df$bfa_env <- sapply(df$bfa_env, function(x) gsub(s_str, "", x))
  }
  
  names(df)[names(df) == source_col] <- "source"
  
  if (all((grepl(remove_str_from_source,df$source) == FALSE))) {
    remove_str_from_source <- paste0("_", dplyr::last(unlist(strsplit(df$source[1], "_"))))
  }
  
  df$source <- sapply(df$source, function(x) gsub(remove_str_from_source, "", x))
  df$bfa_env <- sapply(df$bfa_env, function(x) gsub(remove_str_from_bfa_env, "", x))
  
  return(df)
}


flag_home <- function(df) {
  df$home <- df$bfa_env == df$source
  return(df)
}


# flag_pleiotropy <- function(df, z_cut = 1.96,
#                             source_col = "source",
#                             cluster_col = "cluster",
#                             home_col = "home") {
#   
#   df$fitness_diff <- ((abs(df$s / df$s_se)) > 1.96)
#   df$pleio_pos <- df$fitness_diff == TRUE & df$s > 0
#   df$pleio_neg <- df$fitness_diff == TRUE & df$s < 0
#   
#   suppressWarnings(
#     df %>%
#       dplyr::filter(!!sym(home_col) == FALSE) %>%
#       dplyr::group_by(!!sym(source_col),
#                       !!sym(cluster_col)) %>%
#       dplyr::summarise(s_mu_away = mean(s),
#                        n_bfa_envs = length(unique(bfa_env)),
#                        n_pleio_pos = sum(pleio_pos),
#                        n_pleio_neg = sum(pleio_neg),
#                        n_pleio_tot = sum(fitness_diff)) %>%
#       dplyr::left_join((df %>%
#                           dplyr::filter(!!sym(home_col) == TRUE) %>%
#                           dplyr::group_by(!!sym(source_col),
#                                           !!sym(cluster_col)) %>%
#                           dplyr::summarise(s_mu_home = mean(s))),
#                        by = c(source_col, cluster_col)) ->
#       df_summarized
#   )
#   return(df_summarized)
# }


flag_adapted_at_home <- function(df, source_ref, s_cutoff = 0) {
  
  env_cols <- names(df)[grepl("_s", names(df))]
  
  df_long <-
    df %>%
    tidyr::pivot_longer(cols = all_of(env_cols), names_to = "bfa_env", values_to = "s")
  
  df_long$bfa_env <-  sapply(df_long$bfa_env,
                             function(x) {
                               x%>%
                                 strsplit("\\.") %>%
                                 unlist() %>%
                                 dplyr::first()
                             })
  suppressWarnings(
    df_long_w_flags <- 
      df_long %>%
      dplyr::left_join(source_ref, by = "Full.BC") %>%
      normalize_envs(source_col = "Subpool.Environment",
                     remove_str_from_source = "_2N") %>%
      flag_home()
  )
  
  source_env_check <- with(
    dplyr::filter(df_long_w_flags, home == TRUE),
    table(source, bfa_env)
  )
  
  stopifnot(sum(c(source_env_check[lower.tri(source_env_check)],
                  source_env_check[upper.tri(source_env_check)])) == 0)
  
  df_long_w_flags %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n_home = sum(home)) %>%
    dplyr::filter(n_home == 0) %>% 
    dplyr::select(source) ->
    sources_with_no_homes
  
  df_long_w_flags %>%
    dplyr::filter(home == TRUE, is_adapted == TRUE, s > s_cutoff) ->
    adapted_at_home
  
  df_long_w_flags$adapted_at_home <- FALSE
  df_long_w_flags$adapted_at_home[df_long_w_flags$Full.BC %in% adapted_at_home$Full.BC] <- TRUE
  df_long_w_flags$adapted_at_home[df_long_w_flags$source %in% sources_with_no_homes$source] <- NA
  
  return(dplyr::select(df_long_w_flags, Full.BC, adapted_at_home) %>% unique())
}


load_files_from_arg <- function(arg) {
  arg %>%
    strsplit(" ") %>%
    unlist() ->
    files
  lapply(files, readr::read_csv) ->
    df_list
  
  return(df_list)
}


write_out <- function(df,
                      base_name = NULL,
                      split_on = "_adapted",
                      out_dir = NULL,
                      str_to_append = NULL) {
  
  stopifnot(!is.null(out_dir) & !is.null(str_to_append))
  
  if (!dir.exists(out_dir)) {
    dir.create(file.path(out_dir))
  }
  if (is.null(base_name)) {
    base_name <- "output_"
  } else {
    base_name %>%
      strsplit(split_on) %>%
      unlist() %>%
      dplyr::first() ->
      base_name
  }
  outpath <- file.path(out_dir, paste0(base_name, str_to_append, ".csv"))
  readr::write_csv(df, path = outpath, col_names = T)
}
