#!/usr/bin/env Rscript

# flag_autodiploids.R:
#
# This script ingests hBFA fitness data
# and flags autodiploids based on Which.Subpool identifier passed as argument,
# clusters with distance metric based on traces from known autodips,
# determine appropriate cutoff for assignment to autidip lineages,
# and outputs original input with autodip flag derived from clustering.

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(docopt)))
suppressWarnings(suppressMessages(library(knitr)))
suppressWarnings(suppressMessages(library(kableExtra)))
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path("scripts/src/pleiotropy_functions.R"))

# ------------- #
# function defs #
# ------------- #

# generate distance calculations
distify_autodips <- function(autodips, all_others, cutoff) {

  # calculate distance between each barcode and the set of autodip barcodes
  # sum the total Euclidean distance for each barcode
  # compare this distance for each BC in all_others
  # to the distn created by comparing each autodip to the autodip set.
  # This will generate an empirical distribution,
  # and by applying a cutoff we can determine where to assign all_other BCs
  # to the true autodiploid set.

  # first calculate testing distribution
  autodip_bcs <- row.names(autodips)
  res <- data.frame(NULL)
  
  for (bc in seq_along(autodip_bcs)) {
    ad_tmp   <- autodips[-bc, ]
    focal_bc <- autodips[bc, ]

    # calculate Euclidean distances
    dists <- apply(ad_tmp, 1, function(x) sqrt(sum(x - focal_bc)^2))
    res <- rbind(res,
                 data.frame(Full.BC = autodip_bcs[bc],
                            SQRT_SSED = sum(dists))
    )
  }

  # now calculate distn for all others
  all_other_bcs <- row.names(all_others)
  cat(sprintf("\nnum all_other_bcs:%s", length(all_other_bcs)))
  res_others <- data.frame(NULL)
  
  for (bc in seq_along(all_other_bcs)) {
    focal_bc <- all_others[bc, ]
    
    cat(sprintf("focal other BC:\t%s", focal_bc))
    
    # calculate dist
    dists <- apply(autodips, 1, function(x) sqrt(sum(x - focal_bc)^2))
    cat(sprintf("\nnrow(res_others) pre rbind:\t%s", nrow(res_others)))
    res_others <- rbind(res_others,
                        data.frame(Full.BC = all_other_bcs[bc],
                                   SQRT_SSED = sum(dists)))
    cat(sprintf("\nnrow(res_others) post rbind:\t%s", nrow(res_others)))
  }
  
  cat(sprintf("\nnrow(res) = %s", nrow(res)))
  cat(sprintf("\nnrow(res_others) = %s", nrow(res_others)))
  
  res_tot <- rbind(data.frame(res, autodip = 1),
                   data.frame(res_others, autodip = 0))

  return(res_tot)
}


assign_autodip_status <- function(autodip_dists, cutoff) {

  # define sum SQRT_SSED dist cutoff value
  if (cutoff == "max") {
    d_cutoff <- ceiling(max(autodip_dists$SQRT_SSED[autodip_dists$autodip == 1]))
  } else if (is.numeric(cutoff)) {
    d_cutoff <- ceiling(stats::quantile(autodip_dists$SQRT_SSED[autodip_dists$autodip == 1],
                                 probs = 1 - cutoff)) %>%
    as.numeric()
  }

  # apply cutoff to generate flags; retain barcodes to assign as autodips
  assigned_autodips <- autodip_dists$Full.BC[autodip_dists$autodip == 0 & autodip_dists$SQRT_SSED <= d_cutoff]

  return(list(cutoff = d_cutoff,
              bcs = assigned_autodips)
         )
}


plot_autodip_status <- function(autodip_dists, cutoff, outdir, base_name) {

  n_dp <- sum(autodip_dists$autodip)
  n_other <- nrow(autodip_dists) - n_dp

  autodip_dists <-
    autodip_dists %>%
    dplyr::mutate(ad_flag = ifelse(SQRT_SSED <= cutoff, 1, 0))

  n_ad_flag <- sum(autodip_dists$ad_flag[autodip_dists$autodip == 0])

  # make table
  tab_dat <-
    autodip_dists %>%
    dplyr::filter(autodip == 0) %>%
    dplyr::group_by(Subpool.Environment) %>%
    dplyr::summarise(n_autodip_flags = sum(ad_flag)) %>%
    dplyr::arrange(desc(n_autodip_flags))

  # output table
  tab_dat %>%
    kable(format = "html") %>%
    kable_styling(bootstrap_options = c("striped", "condensed")) %>%
    writeLines(con = file.path(outdir,
                               base_name,
                               "autodip_tally_by_env.html"))

  # make plot
  autodip_dists %>%
    dplyr::mutate(autodip = ifelse(autodip == 1,
                                   "AD+ BCs reference set",
                                   "all other BCs")) %>%
    ggplot() +
    geom_histogram(aes(x = SQRT_SSED, fill = factor(ad_flag)), bins = 50, alpha = 0.8) +
    facet_wrap(~ autodip, scales = 'free_y') +
    scale_fill_manual(values = c("gray40", "darkorange2"), name = "below cutoff") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA),
          plot.title = element_text(size = 9, face = "bold")) +
    geom_vline(xintercept = cutoff, col = "gray60", lty = 2) +
    ggtitle(paste0("Autodip. D distn (",
                   n_dp,
                   " reference BCs, ",
                   n_other,
                   " other BCs)\n",
                   base_name,
                   "\nDist cutoff = ",
                   cutoff,
                   "\nNum. BCs <= cutoff = ",
                   n_ad_flag)) ->
    dist_plot_1

  # output plot
  autodip_dists$Subpool.Environment <- factor(autodip_dists$Subpool.Environment,
                                              levels = tab_dat$Subpool.Environment)

  autodip_dists %>%
    dplyr::filter(autodip == 0) %>%
    ggplot() +
    geom_histogram(aes(x = SQRT_SSED, fill = factor(ad_flag)), bins = 50, alpha = 0.8) +
    facet_wrap(~ Subpool.Environment, nrow = 5, scales = "free_y") +
    scale_fill_manual(values = c("gray40", "darkorange2"), name = "below cutoff") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA)) +
    geom_vline(xintercept = cutoff, col = "gray60", lty = 2) ->
    dist_plot_2

  suppressWarnings(
    ggpubr::ggarrange(plotlist = list(dist_plot_1,dist_plot_2), nrow = 2, align = 'hv',
                    common.legend = T, heights = c(1,3),
                    labels= c("a","b")) %>%
    ggsave(filename = file.path(outdir, base_name, "autodip_tally_plots.png"),
           device = 'png',
           dpi = 300,
           width = 4.5,
           height = 9)
  )
}


# ------------------------- #
# setup command line inputs #
# ------------------------- #

"filter_autodiploids.R

Usage:
    filter_autodiploids.R [--help | --version]
    filter_autodiploids.R [options] <infile> <audodip_tag>

Options:
    -h --help                     Show this screen.
    -v --version                  Show version.
    -o --outdir=<outdir>          Output directory [default: ./]
    -b --base_name=<base_name>    Base name for output file [default: tmp]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -c --cutoff=<pval>            P-value cutoff for outlier trimming [default: 0.05]
    -e --exclude=<env>...         Grep-able string of environments to exclude

Arguments:
    infile                        Input file containing fitness calls for BFA run.
    audodip_tag                   String denoting identifier for autodip lineages
" -> doc

# -------------------- #
# function definitions #
# -------------------- #

run_args_parse <- function(debug_status) {

  if (debug_status == TRUE) {
    arguments <- list()
    arguments$use_iva     <- TRUE
    arguments$infile      <- "data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds.csv"
    arguments$outdir      <- "data/fitness_data/fitness_calls"
    arguments$autodip_tag <- "autodiploids"
    arguments$base_name   <- "hBFA1_cutoff-5"
    arguments$exclude     <- "X48Hr"
    arguments$cutoff      <- 0
    arguments$gens        <- 8
  } else if (debug_status == FALSE) {
    arguments <- docopt::docopt(doc, version = "filter_autodiploids.R v.1.0")
  }
  return(arguments)
}


main <- function(dat, arguments) {

  # filter input dataframe to prepare it for
  cat(sprintf("Preparing fitness file %s...", arguments$infile))
  autodips <-
    dat %>%
    prep_fit_matrix(is_autodip  = TRUE,
                    autodip_tag = arguments$autodip_tag,
                    excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    means_only  = TRUE)

  all_others <-
    dat %>%
    prep_fit_matrix(is_autodip  = FALSE,
                    autodip_tag = arguments$audodip_tag,
                    excludes    = arguments$exclude,
                    iva_s       = arguments$use_iva,
                    gens        = as.double(arguments$gens),
                    means_only  = TRUE)
  
  cat(sprintf("\nnrow(all_others):\t%s", nrow(all_others)))
  
  stopifnot(!is.null(row.names(all_others)), !is.null(row.names(autodips)))
  
  all_others <- all_others[!row.names(all_others) %in% row.names(autodips), ]
  
  cat("Done!\n")
  cat("Flagging autodiploid lineages...")
  autodip_dists <- distify_autodips(autodips   = autodips,
                                    all_others = all_others)

  # now I need to consider the environment for some of these BCs
  # which have on-the-fence distances.
  # those from same Subpool.Env have a higher chance of actually being related
  # if they're not from the ancestral group.

  # need to break this out by subpool.environment
  suppressWarnings(
    autodip_dists %>%
    dplyr::left_join(dat[,c("Full.BC", "Subpool.Environment")],
                     by = "Full.BC") ->
    autodip_dists
  )

  assigned_autodip_output <- assign_autodip_status(autodip_dists, cutoff = arguments$cutoff)

cat("Done!\n")
  cat(sprintf("Generating plots for output in %s...", arguments$outdir))

  plot_autodip_status(autodip_dists,
                      cutoff    = assigned_autodip_output$cutoff,
                      outdir    = arguments$outdir,
                      base_name = arguments$base_name)


  dat$autodip <- FALSE
  dat$autodip[grepl(arguments$autodip_tag, dat$Which.Subpools)] <- TRUE
  dat$autodip[dat$Full.BC %in% assigned_autodip_output$bcs] <- TRUE
  cat("Done!\n")

  return(dat)
}

# ---- #
# main #
# ---- #

debug_status <- FALSE
arguments <- run_args_parse(debug_status)
infile    <- OpenRead(arguments$infile)
dat <- read.table(infile,
                  header = TRUE,
                  sep = ',',
                  stringsAsFactors = F)

cat("\n*************************\n")
cat("* filter_autodiploids.R *\n")
cat("*************************\n\n")

res_out <- main(dat, arguments)

outfile_path <- file.path(arguments$outdir,
                          paste0(arguments$base_name,
                                 "_adapteds_autodips.csv"))

cat(sprintf("Writing output file: %s\n", outfile_path))

write.table(res_out,
            file = outfile_path,
            quote = FALSE,
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

cat("**Script completed successfully!**\n\n")

### END ###
