#! /usr/bin/env Rscript

# call_adapteds.R
#
# this script will take fitness files as input
# and handle a few parameters
# before outputting a table of barcodes and adapted status along with their fitness values.

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(progress)))
suppressWarnings(suppressMessages(library(docopt)))

# ------------------------- #
# setup command line inputs #
# ------------------------- #

'call_adapteds.R

Usage:
    call_adapteds.R [--help | --version]
    call_adapteds.R [options] <infile> <neutral_col>

Options:
    -h --help                     Show this screen.
    -v --version                  Show version.
    -o --outdir=<outdir>          Output directory [default: ./]
    -b --base_name=<name>         Base name for output file [default: tmp]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -c --cutoff=<pval>            P-value cutoff for outlier trimming [default: 0.001]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations
    -r --reps_iter=<sim_reps>     Number of replicate simulations for determining neutral set [default: 500]
    -f --reps_final=<final_reps>  Number of replicate simulations for final neutral probability distribution [default: 1000]

Arguments:
    infile                        Input file containing fitness calls for BFA run.
    neutral_col                   String denoting set of putatively neutral BCs
' -> doc

#print(arguments)

# -------------------- #
# function definitions #
# -------------------- #


run_args_parse <- function(debug_status) {
  
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$use_iva     <- TRUE
    arguments$infile      <- "../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv"
    arguments$outdir     <- "../data/fitness_data/fitness_calls"
    arguments$neutral_col <- 'Ancestor_YPD_2N'
    arguments$name        <- 'cutoff-5'
    arguments$reps_iter   <- 500
    arguments$reps_final  <- 1000
    arguments$exclude     <- 'CLM|FLC4'
    arguments$cutoff      <- 0.01
    arguments$gens        <- 8
  } else if (debug_status == FALSE) {
    arguments <- docopt(doc, version = 'call_adapteds v.1.0')
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

prep_fit_matrix <- function(df, is_neutral = TRUE, means_only = FALSE, neutral_col, excludes, iva_s, gens) {

  # take input counts file and generate list of 2 matrixes:
    # (n x k) matrix of k-variate mean vectors for n barcodes
    # (n x k) k-variate sigma vectors for n barcodes

  # for debugging:
  # excludes <- list(c('CLM'), c('FLC4'), c('Stan'))
  # neutral_col <- 'Ancestor_YPD_2N'
  # df = bfa_dat
  # gens = 8

  avg_type <- ifelse(iva_s == TRUE, 'iva_s', 'ave')

  focal_cols <- c(grep('\\.BC',   names(df), value = T),
                  grep('Subpool', names(df), value = T),
                  grep('\\.R[0-9]', names(df), invert = T, value = T),
                  grep(paste0(avg_type),  names(df), value = T))

  focal_cols <- focal_cols[!grepl(paste0(unlist(excludes), collapse='|'), focal_cols)]

  df <- df[ , names(df) %in% focal_cols]

  # examine list of neutrals if flag present
  if (is_neutral == TRUE) {
    if (!is.null(neutral_col)) {
      df <-
        df %>%
        dplyr::filter(Subpool.Environment == neutral_col)
    } else {
      stop("Please define neutral_col!")
    }
  }

  row.names(df) <- df$Full.BC

  df_mat <-
    df %>%
    dplyr::select(grep('iva_s$|iva_s_err$', names(df), value = T)) %>%
    as.matrix()

  # now I can remove NA, Inf values:
  df_mat[is.infinite(df_mat)] <- NA
  df_mat <- df_mat[complete.cases(df_mat), ]

  # split into means and sigmas:
  if (means_only == FALSE) {

    df_list <- list(means  = df_mat[, grep('iva_s$', colnames(df_mat))] / gens,
                    sigmas = df_mat[, grep('iva_s_err$', colnames(df_mat))] / gens)

    return(df_list)

  } else if (means_only == TRUE) {

    means <- df_mat[, grep('iva_s$', colnames(df_mat))] / gens

    return(means)
  }
}


filter_neutral_outliers <- function(dfs, reps_iter, cutoff, init_cutoff = 0.05, GENS = 8) {

  # for each barcode left out:
  # generate predictive distribution by sampling from vector of means and stderrs

  # # for debugging:
  # dfs <- neut_list
  #
  means  <- dfs[["means"]]
  sigmas <- dfs[["sigmas"]]

  # # grab unique barcodes
  # bcs <- row.names(means) %>% unique()
  #
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

  means  <- means[order(means[,   'starting_dists'], decreasing = TRUE), ]
  sigmas <- sigmas[order(sigmas[, 'starting_dists'], decreasing = TRUE), ]

  # determine quick and dirty the set I'm going to evaluate
  # all with Chisq p-value of < init_cutoff

  init_pvals <- pchisq(means[, 'starting_dists'], df = dim(means)[2] - 1, lower.tail = FALSE)

  # grab barcodes with p-values <= init_cutoff:
  candidate_outliers <- names(init_pvals[init_pvals <= init_cutoff])

  # go through each candidate outlier
  # remove it
  # generate more robust testing distribution
  # flag to remove given p-value

  # NEED TO FINISH THIS ROUTINE:
  # replace with while loop which progressively eliminates barcodes
  # until stopping condition is met
  # for (i in seq_along(candidate_outliers)) {
  #
  #   # remove focal barcode
  #   means2  <- means[-i,  grep('starting_dists', colnames(means), invert = TRUE)]
  #   sigmas2 <- sigmas[-i, grep('starting_dists', colnames(means), invert = TRUE)]
  #   test_bc <- means[i,   grep('starting_dists', colnames(means), invert = TRUE)]
  #
  #   # generate multivariate normal random variables for each vector:
  #   res_mat <- NULL
  #
  #   # wrap this in a function call
  #   for (ii in seq_along(row.names(means2))) {
  #     res_mat <- rbind(res_mat, MASS::mvrnorm(n = R, means2[ii, ], diag(sigmas2[ii, ]^2)))
  #   }
  #
  #   # calculate distance of each barcode in set
  #   D <- mahalanobis(rbind(test_bc,res_mat),
  #                    colMeans(res_mat),
  #                    cov(res_mat))
  #
  #   # generate distance of test barcode
  #   D_test <- mahalanobis(test_bc,
  #                         colMeans(res_mat),
  #                         cov(res_mat))
  #
  #   # calculate p-value of test_bc dist
  #   P_test <- sum(D >= D_test) / length(D)
  # }

  final_neuts <- row.names(means)[!row.names(means) %in% candidate_outliers]
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

main <- function(bfa_dat, arguments) {

  # take bfa_dat
  # run it through name conversion
  cat("Preparing input file...")

  fit_mats <-
    bfa_dat %>%
    prep_fit_matrix(neutral_col = arguments$neutral_col,
                    excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    is_neutral  = TRUE)
  cat("Done!\n")
  cat("Detecting outliers...")

  neutral_set <-
    fit_mats %>%
    filter_neutral_outliers(reps_iter = as.double(arguments$reps_iter),
                            cutoff    = as.double(arguments$cutoff))
  cat("Done!\n")

  # generate testing distribution with refined neutral set
  cat("Generating testing distribution for determining adapteds...")

  neutral_test_distn <-
    fit_mats %>%
    generate_neutral_test_distn(neutrals   = neutral_set,
                                reps_final = as.double(arguments$reps_final))
  cat("Done!\n")

  cat("Determining adapted barcodes...")
  # flag input barcodes as adapted or not based on distance
  adapteds <-
    bfa_dat %>%
    prep_fit_matrix(excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    is_neutral  = FALSE,
                    means_only  = TRUE) %>%
    flag_adapteds(testing_distn = neutral_test_distn,
                  cutoff        = as.double(arguments$cutoff))

  cat("Done!\n")

  # fit_val <- ifelse(arguments$use_iva == TRUE, 'iva_s', 'ave_s')
  # fit_err <- ifelse(arguments$use_iva == TRUE, 'iva_s_err', 'ave_err')
  #
  # fit_cols <- c(grep(paste0(fit_val,'$'), names(bfa_dat), value = T),
  #               grep(paste0(fit_err,'$'), names(bfa_dat), value = T))
  #

  suppressWarnings(
    adapteds_df <-
      bfa_dat %>%
      dplyr::left_join(dplyr::select(adapteds, Full.BC, maha_dist, dist_pval, is_adapted),
                       by = 'Full.BC')
  )

  adapteds_df$neutral_set <- FALSE
  adapteds_df$neutral_set[adapteds_df$Full.BC %in% neutral_set] <- TRUE

  return(adapteds_df)
}

# ---- #
# MAIN #
# ---- #

infile   <- OpenRead(arguments$infile)

# # FOR DEBUGGING
# infile <- OpenRead("../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv")

bfa_dat  <- read.table(infile,
                       header = TRUE,
                       sep = ',',
                       stringsAsFactors = F)

# run
arguments <- run_args_parse(debug_status == TRUE)
res_out <- main(bfa_dat, arguments = arguments)

cat("Writing output file...")

# write
if (!dir.exists(arguments$outdir)) {
  dir.create(arguments$outdir)
}

write.table(res_out,
            file = file.path(arguments$outdir, paste0(arguments$base_name, '_adapteds_',Sys.Date(),'.csv')),
            quote = F,
            sep = ',',
            col.names = T,
            row.names = F)

cat("Done!\n")
cat("**Script completed successfully!**\n\n")

### END ###



#### OBSOLETE BELOW THIS LINE ####

# # grab each environment:
# bfa_envs <-
#   neuts %>%
#   names() %>%
#   grep('Subpool|BC', ., invert = T, value = T) %>%
#   strsplit('\\.') %>%
#   sapply(function(x) dplyr::first(x)) %>%
#   unique()
#
# # run through environments and calculate leave-one-out z-score in each environment
# res_by_env <- lapply(bfa_envs, function(x) detect_outliers(env = x, df = neuts))
#
# # bring scores per barcode together and generate plot of distributions
# res_by_env_long <- do.call(rbind, res_by_env)
#
# # generate sum of squares to see if it has decent testing properties
# res_by_env_long %>%
#   dplyr::group_by(Full.BC, bfa_env) %>%
#   dplyr::summarise_if(is.numeric, sum) ->
#   res_by_env_sum
#
# # the squared z-score for each environment should be distributed as a Chi-squared with 1 df.
# # we can therefore calculate the joint likelihood across environments and apply a cutoff:
#
#
# res_by_env_long %>%
#   dplyr::group_by(Full.BC, bfa_env) %>%
#   dplyr::summarise(pval = pchisq(q0.500, df = 1, lower.tail = FALSE)) ->
#   res_by_env_p
#
# res_by_env_p_prod <-
#   res_by_env_p %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(Full.BC) %>%
#   dplyr::summarise(pprod = prod(pval, na.rm = T)) %>%
#   dplyr::full_join(
#     tidyr::spread(res_by_env_p, key = 'bfa_env', value = 'pval'))
#
# ggplot(res_by_env_p_prod, aes(x = pprod)) +
#   stat_ecdf() #geom_histogram(bins = 100)
#
#
#
# #### obsolete below this line ####
#
#
# ggplot(res_by_env_sum, aes(x = q0.500)) + #geom_histogram(bins = 50) +
#   geom_density() +
#   #facet_wrap(~ bfa_env) +
#   scale_x_continuous(limits = c(0,3))
#
# # try plotting
# ggplot(res_by_env_long, aes(x = bfa_env, y = q0.500, ymin = q0.025, ymax = q0.975, group = Full.BC)) +
#   geom_linerange(alpha = 0.1) +
#   geom_line(alpha = 0.3) +
#   geom_point(alpha = 0.5) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
#   scale_y_log10()
#
# # do sum across environment and see what it looks like
# res_by_env_sum <-
#   res_by_env_long %>%
#   dplyr::group_by(Full.BC) %>%
#   dplyr::summarise_if(is.numeric, sum)
#
# # need to define cutoff here somehow.
# # In theory, this will be the same process as for defining adapted status, but with the Neutrals.
# # calculate p-value under the Chisquared with df = 1 for each z_i
# # inspect distribution of log p-values across all environments; find the product across environments
# # Basically calculating the likelihood of the data under the hypothesis that mean(X) = 0.
#
# ggplot(res_by_env_sum, aes(x = q0.500)) + geom_density(bins = 1000) +
#   scale_x_continuous(limits = c(0,50))
#
# # now try to implement:
# # run on each environment; decide on cutoff for excluding outliers
# res_test <- dplyr::arrange(res_test, q0.500)
# res_test$Full.BC <- factor(res_test$Full.BC, levels = paste0(res_test$Full.BC))
#
# ggplot(res_test, aes(x = Full.BC, y = q0.500)) +
#   geom_point() +
#   geom_linerange(aes(ymin = q0.025, ymax = q0.975)) +
#   theme_minimal()
