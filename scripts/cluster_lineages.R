#! /usr/bin/env Rscript

# cluster_lineages.R
#
# this script will take input fitness files 
# and generate t-SNE based clustering 
# with exclusions (hard-coded for now)

# ------------------------- #
# header                    #
# ------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(progress)))
suppressWarnings(suppressMessages(library(Rtsne)))
suppressWarnings(suppressMessages(library(docopt)))
source(file.path('./src/adapted_functions.R'))

# ------------------------- #
# setup command line inputs #
# ------------------------- #

'cluster_lineages.R

Usage:
    cluster_lineages.R [--help | --version]
    cluster_lineages.R [options] INFILE ...

Options:
    -h --help                     Show this screen.
    -v --version                  Show version.
    -o --outdir=<outdir>          Output directory [default: ./]
    -u --use_iva                  Flag to determine whether to use inverse variance weighted avg or arithmentic avg [default: TRUE]
    -g --gens=<gens>              Number of generations per cycle (used to divide input fitness estimates) [default: 8]
    -e --exclude=<env>...         Space-separated list of environments to exclude from neutral set calculations
    
Arguments:
    INFILE                        Input file(s) containing fitness calls for BFA run.
' -> doc

# -------------------- #
# function definitions #
# -------------------- #

run_args_parse <- function(debug_status) {
  
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$use_iva     <- TRUE
    arguments$infile      <- c("../data/fitness_data/fitness_calls/dBFA2_cutoff-5_adapteds_2020-03-03.csv","../data/fitness_data/fitness_calls/hBFA1_cutoff-5_adapteds_autodips.csv")
    arguments$outdir      <- "../data/fitness_data/fitness_calls/clusters"
    arguments$exclude     <- 'CLM|FLC4|Stan|X48Hr|X02M'
    arguments$gens        <- 8
  } else if (debug_status == FALSE) {
    arguments <- docopt(doc, version = 'cluster_lineages.R v.1.0')
  }
  return(arguments)
}

filter_to_focal_bcs <- function(x,
                           adapt_col   = 'is_adapted', 
                           neutral_col = 'neutral_set',
                           autodip_col = 'autodip',
                           retain_neutrals = TRUE) {
  
  # goal: retain neutral sets.
  if (autodip_col %in% names(x)) {
    x2 <- x$Full.BC[(x[, adapt_col] == TRUE) & (x[, autodip_col] == FALSE)]
    x2  <- c(x2, x$Full.BC[x[, neutral_col] == TRUE])
  } else {
    x2 <- x$Full.BC[x[, adapt_col] == TRUE]
  }
  
  if (retain_neutrals == TRUE) {
    x2  <- c(x2, x$Full.BC[x[, neutral_col] == TRUE])
  }
  return(x2)
}


check_envs <- function(x) {
  
  # make sure envs are shared between all files
  the_names <- sapply(x, function(x) sort(colnames(x)), simplify = FALSE)
  
  # test out
  intersect_mat <- sapply(the_names,
                          function(x) sapply(the_names, 
                                             function(y) length(intersect(x,y))))
  # zero out diag elements
  intersect_mat <- intersect_mat * (1 - diag(dim(intersect_mat)[1]))
  
  # sum off-diagonal.
  off_diag_sum = sum(intersect_mat, na.rm = T)
  
  # check for correct sum
  target_sum <- max(sapply(the_names, length)) * length(the_names)
  
  if (off_diag_sum == target_sum) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

do_cluster <- function(dat, meta) {
  
  set.seed(12345)
  
  tsne.norm = Rtsne(dat, pca = FALSE)
  
  return(tsne.norm)
}

grab_metadata <- function(x,
                          bcs, 
                          target_names = c('Subpool.Environment','Which.Subpools', 'neutral_set')) {
  
  x2 <- x[x$Full.BC %in% bcs, names(x) %in% c('Full.BC', target_names)]
  
  # parse to source_env and ploidy
  x2$ploidy     <- sapply(x2$Subpool.Environment, function(x) ifelse(grepl('2N',x),'2N','1N'))
  x2$source_env <- sapply(x2$Subpool.Environment,
                          function(s) {
                            tmp <- unlist(strsplit(s, '_'))
                            paste0(tmp[1:(length(tmp) - 1)], collapse = '_')
                          })
  
  x2$source_env <- sapply(x2$source_env, function(x) ifelse(grepl('YPD',x),'YPD',x))
  
  return(x2)
}

main <- function(arguments) {
  
  # check + make output dir
  if (!dir.exists(arguments$outdir)) {
    dir.create(file.path(arguments$outdir))
  }
  
  # check + grab input files
  if (sum(!sapply(arguments$infile, file.exists)) > 0 ) {
    stop("One or more infile does not exist. Please check infile path and try again.", call. = FALSE)
  }
  
  # read infiles
  #infiles <- sapply(arguments$infile, function(x) OpenRead(file.path(x)))
  infiles <- sapply(arguments$infile, function(x) read.table(file.path(x),
                                                                   sep = ',',
                                                                   header = T,
                                                                   stringsAsFactors = F))
  
  # filter out neutrals and autodips
  
  
  # prep the matrixes:
  fit_prepped <- lapply(infiles,
                        prep_fit_matrix,
                        means_only = TRUE,
                        excludes = arguments$exclude,
                        iva_s = arguments$use_iva,
                        gens = arguments$gens)
  
  # grab barcodes to retain after filtering:
  bcs_to_retain <- lapply(infiles, filter_to_focal_bcs, retain_neutrals = TRUE)
  bcs_to_retain <- do.call(c, bcs_to_retain)
  
  # check for only share envs among infiles
  if (!check_envs(fit_prepped)) {
    stop("Environments are not shared between all input files. Exiting.", call. = FALSE)
  }
  
  # combine fit_prepped; initiate clustering
  dat <- do.call(rbind, fit_prepped)
  dat <- dat[row.names(dat) %in% c(bcs_to_retain), ]
  
  cdat <- do_cluster(dat)
  info.norm <- data.frame(tsne1 = cdat$Y[, 1],
                          tsne2 = cdat$Y[, 2])
  
  # grab meta-data from original DFs
  meta <- lapply(infiles, grab_metadata, bcs = row.names(dat)) %>%
    do.call(rbind, .)
  
  # join; send to plot
  dat <-
    dat %>%
    data.frame()
  
  dat$Full.BC <- row.names(dat)
  
  dat_full <-
    dat %>%
    dplyr::left_join(meta, by = 'Full.BC') %>%
    dplyr::select(Full.BC, everything()) %>%
    dplyr::mutate(tsne1 = cdat$Y[, 1],
                  tsne2 = cdat$Y[, 2])
  
  # plot t-sne; need to add meta-data for clustering
  ggplot(dat_full, aes(x = tsne1, y = tsne2, col = source_env)) + 
    geom_point(alpha = 0.3) + theme_bw() +
    facet_wrap(~ ploidy)
    
  
  # need to find various means of clustering; bring in mutation data; 
  # look for mutual information.
}

debug_status <- TRUE
arguments <- run_args_parse(debug_status)
#print(arguments)
main(arguments)

