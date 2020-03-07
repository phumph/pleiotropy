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
source(file.path('./src/adapted_functions.R'))

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
    -b --base_name=<base_name>    Base name for output file [default: tmp]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -c --cutoff=<pval>            P-value cutoff for outlier trimming [default: 0.05]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations
    -r --reps_iter=<sim_reps>     Number of replicate simulations for determining neutral set [default: 1000]
    -f --reps_final=<final_reps>  Number of replicate simulations for final neutral probability distribution [default: 1000]

Arguments:
    infile                        Input file containing fitness calls for BFA run.
    neutral_col                   String denoting set of putatively neutral BCs
' -> doc

# -------------------- #
# function definitions #
# -------------------- #

run_args_parse <- function(debug_status) {

  if (debug_status == TRUE) {
    arguments <- list()
    arguments$use_iva     <- TRUE
    arguments$infile      <- "../data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv"
    arguments$outdir      <- "../data/fitness_data/fitness_calls"
    arguments$neutral_col <- 'Ancestor_YPD_2N'
    arguments$base_name   <- 'dBFA2_cutoff-5'
    arguments$reps_iter   <- 1000
    arguments$reps_final  <- 1000
    arguments$exclude     <- 'CLM|FLC4|Stan'
    arguments$cutoff      <- 0.05
    arguments$gens        <- 8
  } else if (debug_status == FALSE) {
    arguments <- docopt(doc, version = 'call_adapteds v.1.0')
  }
  return(arguments)
}

main <- function(bfa_dat, arguments) {

  set.seed(12345)

  # take bfa_dat
  # run it through name conversion
  cat("Preparing input file...")

  fit_mats <-
    bfa_dat %>%
    prep_fit_matrix(neutral_col = arguments$neutral_col,
                    excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    is_neutral  = TRUE,
                    do_filter   = TRUE,
                    means_only  = FALSE)
  cat("Done!\n")
  cat("Detecting outliers in neutral barcode set...")

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
  # write to file:
  if (!dir.exists(file.path(arguments$outdir, arguments$base_name))) {
    dir.create(file.path(arguments$outdir, arguments$base_name))
  }

  saveRDS(neutral_test_distn,
          file = file.path(arguments$outdir, arguments$base_name, 'neutral_test_distn.Rds'))

  # save plot of distribution
  data.frame(distances = neutral_test_distn$distances) %>%
    ggplot() +
    geom_histogram(aes(x = distances), bins = 100) +
    theme_minimal() +
    theme(axis.line = element_line(),
          plot.title = element_text(size = 9, face = 'bold')) +
    geom_vline(xintercept = quantile(neutral_test_distn$distances, probs = 1 - as.double(arguments$cutoff)), col = 'darkorange') +
    ggtitle(paste0("Neutral D distn (", length(neutral_set)," BCs)\n", arguments$base_name, "\ncutoff = ", arguments$cutoff)) ->
    neutral_dist_plot

  ggsave(neutral_dist_plot,
         filename = file.path(arguments$outdir, arguments$base_name, 'neutral_test_distn.png'),
         device = 'png',
         dpi = 300,
         width = 5,
         height = 3.5,
         units = 'in')

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

# run
cat("*******************\n")
cat("* call_adapteds.R *\n")
cat("*******************\n\n")

debug_status <- FALSE
arguments <- run_args_parse(debug_status)
infile    <- OpenRead(arguments$infile)
dat <- read.table(infile,
                 header = TRUE,
                 sep = ',',
                 stringsAsFactors = F)

res_out <- main(dat, arguments)

# write
if (!dir.exists(arguments$outdir)) {
  dir.create(arguments$outdir)
}

outfile_path <- file.path(arguments$outdir, paste0(arguments$base_name, '_adapteds.csv'))

cat(sprintf("Writing output file: %s\n", outfile_path))

write.table(res_out,
            file = outfile_path,
            quote = F,
            sep = ',',
            col.names = T,
            row.names = F)

cat("**Script completed successfully!**\n\n")

### END ###
