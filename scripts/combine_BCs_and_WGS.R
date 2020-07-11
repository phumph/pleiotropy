#!/usr/bin/env Rscript

# script to take BC file and add BC info to WGS mutation calls.

# inputs (hard-coded):
  # 1. Barcodes file
  # 2. Mutation calls file
  # 3. File containing all div--env BC pairs from BFAs
  # 4. Flag for rules on how to call majority BCs

### >>>>>>
### Header
### >>>>>>

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))

f_bfa_bcs  <- file.path('data/mutation_data/barcode_extraction/PLT_all_BFA_bcs.csv')
f_wgs_bcs  <- file.path('data/mutation_data/barcode_extraction/bcs_extracted_compiled.csv')
f_variants <- file.path('data/mutation_data/master_mutation_calls.txt')
f_subpools_dbfa2 <- file.path('data/fitness_data/bfa_bc_counts/dBFA2_counts_with_env_info.csv')
f_subpools_hbfa1 <- file.path('data/fitness_data/bfa_bc_counts/hBFA1_counts_with_env_info.csv')
f_subpools_hbfa2 <- file.path('data/fitness_data/bfa_bc_counts/hBFA2_counts_with_env_info.csv')

bfa_bcs  <- read.table(f_bfa_bcs,
                       sep = ',', header = T, stringsAsFactors = F, na.strings = c('',' '))
wgs_bcs  <- read.table(f_wgs_bcs,
                       sep = ',', header = T, stringsAsFactors = F, na.strings = c('',' '))
variants <- read.table(f_variants,
                       sep = '\t', header = T, stringsAsFactors = F, na.strings = c('',' '),
                       allowEscapes = FALSE, skipNul = TRUE)
subpools_dbfa2 <- read.table(f_subpools_dbfa2,
                             sep = ',', header = T, stringsAsFactors = F, na.strings = c('',' '))
subpools_hbfa1 <- read.table(f_subpools_hbfa1,
                             sep = ',', header = T, stringsAsFactors = F, na.strings = c('',' '))
subpools_hbfa2 <- read.table(f_subpools_hbfa2,
                             sep = ',', header = T, stringsAsFactors = F, na.strings = c('',' '))

### >>>>>>><<<<<<<<>>>>>>>>><<<<<<<<<>>>>>>
### Reconciling strain and file identifiers
### >>>>>>><<<<<<<<>>>>>>>>><<<<<<<<<>>>>>>

# define proportion cutoff for calling BCs
CUTOFF <- 0.8

# helper function to take each row and find its matches
add_bc_to_wgs <- function(x, ...) {
  # find matching row
  the_matches <- grep(x, wgs_bcs$File)
  # BC grabbing routine
  if (length(the_matches) > 0) {
    res_full <- data.frame()
    # add for loop to cycle through all of the matches
    for(b in seq_along(the_matches)) {
      the_match <- the_matches[b]
      if (!grepl(';', wgs_bcs$dbcs[the_match]) & !grepl(';', wgs_bcs$ebcs[the_match])) {
        Diverse.BC           <- wgs_bcs$dbcs[the_match]
        Diverse.BC_count     <- as.numeric(wgs_bcs$dbc_counts[the_match])
        Environment.BC       <- wgs_bcs$ebcs[the_match]
        Environment.BC_count <- as.numeric(wgs_bcs$ebc_counts[the_match])
        dbc_prop <- ifelse(!is.na(wgs_bcs$dbcs[the_match]), 1, NA)
        ebc_prop <- ifelse(!is.na(wgs_bcs$ebcs[the_match]), 1, NA)
      } else if (grepl(';', wgs_bcs$dbcs[the_match]) | grepl(';', wgs_bcs$ebcs[the_match])) { # case where more than one count in either BC type
        cat(sprintf("File with multiple barcodes: %s\n", x))
        # grab most common dbc if frequency is > 20% of total counts
        dbc_counts <- strsplit(wgs_bcs$dbc_counts[the_match], ';') %>% unlist() %>% as.numeric()
        dbc_prop <- dbc_counts[1] / sum(dbc_counts)
        if (!is.na(dbc_prop) & dbc_prop >= CUTOFF) {
          Diverse.BC <- strsplit(wgs_bcs$dbcs[the_match], ';') %>% unlist() %>% dplyr::first()
          Diverse.BC_count <- dbc_counts[1]
        } else {
          # do something else if doesn't meet CUTOFF
          # for now, return INS for insufficient coverage
          Diverse.BC <- 'INS'
          Diverse.BC_count <- dbc_counts[1]
        }
        ebc_counts <- strsplit(wgs_bcs$ebc_counts[the_match], ';') %>% unlist() %>% as.numeric()
        ebc_prop <- ebc_counts[1] / sum(ebc_counts)
        if (!is.na(ebc_prop) & ebc_prop >= CUTOFF) {
          Environment.BC <- strsplit(wgs_bcs$ebcs[the_match], ';') %>% unlist() %>% dplyr::first()
          Environment.BC_count <- ebc_counts[1]
        } else {
          # do something else if doesn't meet CUTOFF
          # for now, return INS for insufficient coverage
          Environment.BC <- 'INS'
          Environment.BC_count <- ebc_counts[1]
        }
      }
      res <- data.frame(Strain = x,
                        BC_Strain_ID = wgs_bcs$File[the_match],
                        Diverse.BC,
                        Environment.BC,
                        Diverse.BC_count,
                        Environment.BC_count,
                        dbc_prop = round(dbc_prop,2),
                        ebc_prop = round(ebc_prop,2),
                        stringsAsFactors = FALSE)
      res_full <- dplyr::bind_rows(res_full, res)
    }
  } else {
    # no match:
    res_full <- data.frame(Strain = x,
                           BC_Strain_ID = NA,
                           Diverse.BC = NA,
                           Environment.BC = NA,
                           Diverse.BC_count = NA,
                           Environment.BC_count = NA,
                           dbc_prop = NA,
                           ebc_prop = NA)
  }
  return(res_full)
}

sapply(variants$Strain, function(x) sum(grepl(tolower(paste0(x)), tolower(wgs_bcs$File)))) -> t1
variants[variants$Strain %in% names(t1[t1==1]), ] -> vars_1x
variants[variants$Strain %in% names(t1[t1==2]), ] -> vars_2x

# build join key for vars_1x:
l1 <- as.list(unique(vars_1x$Strain))

# safe to ignore warnings
lapply(l1, function(x) add_bc_to_wgs(x)) %>%
  do.call(rbind, .) ->
  vars_1x_BCs

# build join key for vars_2x
# safe to ignore warnings
l2 <- as.list(unique(vars_2x$Strain))

lapply(l2, function(x) add_bc_to_wgs(x)) %>%
  do.call(rbind, .) ->
  vars_2x_BCs

###
### Merging barcodes with variant calls
###

# need to filter out barcodes with ambiguous calls
vars_1x_filt <-
  vars_1x_BCs %>%
  dplyr::filter(!is.na(Diverse.BC),
                !is.na(Environment.BC),
                dbc_prop >= CUTOFF,
                ebc_prop >= CUTOFF) #%>% nrow() # 586/876 [67%]

vars_2x_filt <-
  vars_2x_BCs %>%
  dplyr::filter(!is.na(Diverse.BC),
                !is.na(Environment.BC),
                dbc_prop >= CUTOFF,
                ebc_prop >= CUTOFF) #%>% nrow() # 59/156 [38%]

# make sure BC calls for redundant strain_ids match:
table(vars_2x_filt$Strain)[table(vars_2x_filt$Strain)>1] %>% names() -> dup_strains

###
### Filtering based on sub-pool environment calls
###

bc_cols <- c('Full.BC',
  'Diverse.BC',
  'Environment.BC',
  'Total.Counts',
  'Subpool.Environment',
  'Which.Subpools',
  'Putative.Environment')

dplyr::bind_rows(
  dplyr::select(subpools_dbfa2, bc_cols),
  dplyr::select(subpools_hbfa1, bc_cols),
  dplyr::select(subpools_hbfa2, bc_cols)
) -> bcs_envs

all_vars <-
  dplyr::bind_rows(vars_1x_filt, vars_2x_filt) %>%
  dplyr::mutate(is_dup = Strain %in% dup_strains) %>%
  dplyr::left_join(bcs_envs, by = c('Diverse.BC','Environment.BC'))

# extract environment information from strain ID
all_vars$wgs_env <- sapply(all_vars$BC_Strain_ID,
                           function(x) {
                             if (grepl('Dip_clone', x)) {
                               strsplit(x, '_') %>%
                                 unlist() -> x2
                               x2[6]
                             } else {
                               strsplit(x, '_') %>%
                                 unlist() %>%
                                 first()
                             }
                           }
)

all_vars$wgs_env[grep('[0-9]', all_vars$wgs_env)] <- gsub('[a-zA-Z]','',all_vars$wgs_env[grep('[0-9]', all_vars$wgs_env)])

# table(all_vars$wgs_env,all_vars$Subpool.Environment)

# annotate numeric environments:
all_vars$wgs_env[all_vars$wgs_env ==  4] <- '37C'
all_vars$wgs_env[all_vars$wgs_env ==  6] <- 'pH7_3'
all_vars$wgs_env[all_vars$wgs_env == 10] <- '21C'
all_vars$wgs_env[all_vars$wgs_env == 11] <- '02M_NaCl'

###
### flag BCs based on environment match
###

all_vars$env_match <- mapply(function(x,y) grepl(x,y), x = all_vars$wgs_env, y = all_vars$Subpool.Environment)

# exclude non-matching
all_vars_matched <-
  all_vars %>%
  dplyr::filter(env_match == TRUE)

###
### bring in barcode information to the variant data
###

variants <-
  variants %>%
  dplyr::left_join(
    dplyr::select(all_vars_matched, Strain, Full.BC, Diverse.BC, Environment.BC),
    by = "Strain")

# filter out those without matches
variants_bc <-
  variants %>%
  dplyr::filter(!is.na(Full.BC))

###
### Export data
###

write.table(variants_bc,
            file = file.path("data/mutation_data/mutations_by_bc.csv"),
            sep = ",",
            quote = T,
            row.names = F,
            col.names = T)
