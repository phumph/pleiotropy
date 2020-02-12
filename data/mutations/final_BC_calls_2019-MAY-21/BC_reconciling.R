# barcode reconciling script
# goal: count up instances of correspondnce and disagreement by source and ploidy; output unique counts
# last updated: 2019-MAY-21 by PTH

library(dplyr)
library(ggplot2)
setwd("~/Dropbox/PLT/analysis_PLT/PLT/data/mutations/final_BC_calls_2019-MAY-21")

# load data
bcs <- read.table("2019-MAY-21_reconciling.txt",T,"\t")

# make call to flag if absent in LH but present in LY
bcs$YL_not_LH <- (bcs$LH_Diverse.BC %in% c("ABSENT","None")) & (!bcs$YL_Diverse.BC %in% c('None'))
#table(bcs$YL_not_LH)

# make call to flag if absent in YL but present in LH
bcs$LH_not_YL <- (!bcs$LH_Diverse.BC %in% c("ABSENT",'None')) & (bcs$YL_Diverse.BC %in% c('None'))
#table(bcs$LH_not_YL)

# make call to flag if barcode in both and they are the same
bcs$DBC_match <- (as.vector(bcs$LH_Diverse.BC) == as.vector(bcs$YL_Diverse.BC)) & (bcs$LH_Diverse.BC != "None")
#table(bcs$DBC_match)

bcs$DBC_conflict <- (as.vector(bcs$LH_Diverse.BC) != as.vector(bcs$YL_Diverse.BC)) & ((!bcs$LH_Diverse.BC %in% c('ABSENT','None')) & (!bcs$YL_Diverse.BC %in% c('None')))
#table(bcs$DBC_conflict)
#bcs[bcs$DBC_conflict==TRUE,grep('Diverse.BC',names(bcs))]

# make call to flag if absent in both:
bcs$no_BC <- (bcs$LH_Diverse.BC %in% c("ABSENT",'None')) & (bcs$YL_Diverse.BC %in% c('None'))
table(bcs$no_BC)

# exclude rows with multiple BC calls from Milo's script
bcs$multiple_DBCs <- FALSE
bcs$multiple_DBCs[grep(';',bcs$YL_Diverse.BC)] <- TRUE
#table(bcs$multiple_DBCs)

# write intermediate file for further reconciling elsewhere:
write.csv(bcs, file = './bcs_transformed_2019-JUN-18.csv',row.names = F)

# exclude rows with multiple barcodes:
bcs2 <- dplyr::filter(bcs,multiple_DBCs==FALSE)

# need to check that there is a flag in every row:
bcs2$flag_sum <- apply(bcs2[,names(bcs2) %in% c('YL_not_LH','LH_not_YL','DBC_match','no_BC','DBC_conflict')], 1, sum)
#table(bcs2$flag_sum)

# find total number of WGS per barcode call per ploidy and source:
WGS_summary <- dplyr::group_by(bcs2, SOURCE, PLOIDY) %>%
  summarise(n_YL_not_LH = sum(YL_not_LH, na.rm = T),
            n_HL_not_YL = sum(LH_not_YL, na.rm = T),
            n_DBC_match = sum(DBC_match, na.rm = T),
            n_DBC_conflict = sum(DBC_conflict, na.rm = T),
            n_no_BC = sum(no_BC, na.rm = T),
            tot_with_DBC = sum(n_YL_not_LH,n_HL_not_YL,n_DBC_match),
            tot_WGS = length(FULL_WGS_SAMPLE_ID),
            prop_WGS = round(tot_with_DBC/tot_WGS,2))

# export table:
write.table(WGS_summary, file = "WGS_summary.txt", row.names = F, quote = F, sep = '\t')

# generate table where we consider unique barcodes.
# construct Diverse.BC column and then tabulate how many WGS samples per barcode (>=2)
# and then how many unique genotypes we have WGS data per env and ploidy

# bring in relevant rows to compile consensus DBCs
bcs2$Diverse.BC <- "no_consensus"
bcs2$Diverse.BC[bcs2$LH_not_YL==TRUE] <- paste0(bcs2$LH_Diverse.BC[bcs2$LH_not_YL==TRUE])
bcs2$Diverse.BC[bcs2$YL_not_LH==TRUE] <- paste0(bcs2$YL_Diverse.BC[bcs2$YL_not_LH==TRUE])
bcs2$Diverse.BC[bcs2$DBC_match==TRUE] <- paste0(bcs2$YL_Diverse.BC[bcs2$DBC_match==TRUE])

# now, grab all unique Diverse.BCs
bcs3 <- dplyr::filter(bcs2, Diverse.BC != 'no_consensus')

# generate table of all uniques
unique_bcs <- dplyr::select(bcs3, Diverse.BC, SOURCE, PLOIDY) %>%
  unique(.) %>% group_by(SOURCE, PLOIDY) %>%
  summarise(n_unique_bcs = length(Diverse.BC)) %>%
  arrange()

# export this table; load into presentation
write.table(unique_bcs, file = "unique_bcs_table.txt", row.names = F, quote = F, sep = '\t')

bc_counts <- dplyr::group_by(bcs3, Diverse.BC, SOURCE, PLOIDY) %>%
  summarise(n_wgs_samples = length(Diverse.BC)) %>% 
  #filter(n_wgs_samples > 1) %>%
  arrange(desc(n_wgs_samples))

# figure out easy way to plot these data by source:
wgs_sample_coverage_p1 <- ggplot(dplyr::filter(bc_counts, n_wgs_samples > 1), aes(x = n_wgs_samples, fill = PLOIDY)) + 
  geom_histogram(bins = 50) + facet_wrap(~ SOURCE, scales='free') +
  scale_fill_manual(values = c('darkorange2','dodgerblue')) +
  theme_bw() + theme(legend.position = "top")

ggsave(wgs_sample_coverage_p1, file = "wgs_sample_coverage_p1.pdf", width = 5, height = 3.5)

# export this table
write.table(bc_counts, file = "n_wgs_samples_count.txt", row.names = F, quote = F, sep = '\t')

# possible next exercise: take all samples with same consensus barcode and generate list of .fastq files to collapse
# to re-run variant calling to increase coverage.