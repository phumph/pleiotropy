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
suppressWarnings(suppressMessages(library(docopt)))

source(file.path('scripts/src/pleiotropy_functions.R'))

# -------------------- #
# function definitions #
# -------------------- #

main <- function(arguments) {
  
  dat <- read.table(arguments$infile,
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = F)
  
  set.seed(12345)
  
  # take bfa_dat
  # run it through name conversion
  cat("Preparing input file...")
  
  fit_mats <-
    dat %>%
    prep_fit_matrix(neutral_col = arguments$neutral_col,
                    excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    is_neutral  = TRUE,
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
  
  adapteds <-
    dat %>%
    prep_fit_matrix(excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    is_neutral  = FALSE,
                    means_only  = TRUE) %>%
    flag_adapteds(testing_distn = neutral_test_distn,
                  cutoff        = as.double(arguments$cutoff))
  
  suppressWarnings(
    dat %>%
      prep_fit_matrix(iva_s       = arguments$use_iva,
                      gens        = as.double(arguments$gens),
                      is_neutral  = FALSE,
                      means_only  = TRUE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::rename("Full.BC" = "rowname") %>%
      dplyr::left_join(
        dplyr::select(adapteds, Full.BC, maha_dist, dist_pval, is_adapted),
        by = "Full.BC") %>%
      dplyr::filter(is_adapted == TRUE) ->
      dat_for_home
    )
  
  source_ref <- dplyr::select(dat, Full.BC, Subpool.Environment)
  adapted_at_home <- flag_adapted_at_home(dat_for_home, source_ref = source_ref)
  
  adapted_at_home %>%
    dplyr::filter(adapted_at_home == FALSE) ->
    non_adapted_BCs
  
  adapteds$is_adapted[adapteds$Full.BC %in% non_adapted_BCs$Full.BC] <- FALSE
  
  cat("Done!\n")
  
  suppressWarnings(
    adapteds_df <-
      dat %>%
      dplyr::left_join(dplyr::select(adapteds,
                                     Full.BC, maha_dist, dist_pval, is_adapted),
                       by = "Full.BC")
  )
  
  adapteds_df$neutral_set <- FALSE
  adapteds_df$neutral_set[adapteds_df$Full.BC %in% neutral_set] <- TRUE
  
  if (!dir.exists(arguments$outdir)) {
    dir.create(arguments$outdir)
  }
  
  outfile_path <- file.path(arguments$outdir,
                            paste0(arguments$base_name,
                                   "_adapteds.csv"))
  
  cat(sprintf("Writing output file: %s\n", outfile_path))
  
  write.table(adapteds_df,
              file = outfile_path,
              quote = F,
              sep = ",",
              col.names = T,
              row.names = F)
}

# ==== #
# main #
# ==== #

"call_adapteds.R

Usage:
    call_adapteds.R [--help]
    call_adapteds.R [options] <infile> <neutral_col>

Options:
    -h --help                     Show this screen.
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
" -> doc

args <- list(
  use_iva     = TRUE,
  infile      = "data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv",
  outdir      = "data/fitness_data/fitness_calls",
  neutral_col = "Ancestor_YPD_2N",
  base_name   = "dBFA2_cutoff-5",
  reps_iter   = 5000,
  reps_final  = 5000,
  exclude     = "CLM|FLC4|Stan",
  cutoff      = 0.05,
  gens        = 8
)

debug_status <- FALSE

cat("*******************\n")
cat("* call_adapteds.R *\n")
cat("*******************\n\n")

arguments <- run_args_parse(args, debug_status)
main(arguments)

cat("**Script completed successfully!**\n\n")
