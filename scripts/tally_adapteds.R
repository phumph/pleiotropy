# tally_adapteds.R

library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
library(ggpubr)

# script takes two inputs as the outputs from:
  # cluster_adapteds.R
  # combine_BCs_and_WGS.R

# and produces plots and tables of counts of adapted barcodes by env with WGS
# this script also performs clustering in order to display WGS counts by cluster

#### DEBUGGING PANEL
arguments <- list()
arguments$use_iva <- TRUE
arguments$adapteds <- '../data/fitness_data/fitness_calls/cutoff-5_adapteds_2019-11-25.csv'
arguments$variants <- '../data/mutation_data/mutations_by_bc.csv'
arguments$exclude  <- 'Stan'
arguments$gens     <- 8

#### MAIN ####

# bring in adapteds
adapteds <- read.table(file = file.path(arguments$adapteds),
                       header = T, 
                       sep = ',',
                       stringsAsFactors = F)

# subset columns to focal_cols
fit_val <- ifelse(arguments$use_iva == TRUE, 'iva_s', 'ave_s')
focal_cols <- c('Full.BC','Diverse.BC','Environment.BC','Subpool.Environment','Which.Subpools',
                grep(fit_val, names(adapteds), value = T),
                'maha_dist','dist_pval','is_adapted','neutral_set')
adapteds <- adapteds[, names(adapteds) %in% focal_cols]


# bring in variants
variants <- read.table(file = file.path(arguments$variants),
                       header = T, 
                       sep = ',',
                       stringsAsFactors = F)


# add flag whether we have WGS data
adapteds$wgs <- FALSE
adapteds$wgs[adapteds$Full.BC %in% variants$Full.BC] <- TRUE

#### PLOTS ####
ttt <- data.frame(table(adapteds$is_adapted,
                        adapteds$Subpool.Environment)) %>% 
  dplyr::arrange(desc(Var1),desc(Freq))

# scrub levels that don't make sense; re-order:
ttt <- ttt[!ttt$Var2 %in% c('none','not_read'),]
ttt$Var2 <- factor(ttt$Var2, levels = paste0(unique(ttt$Var2)))

# supply text annotations:
ttt2 <- reshape2::dcast(ttt, Var2 ~ Var1, value.var = 'Freq')
names(ttt2) <- c('source','non_adapted','adapted')
ttt2 <- dplyr::mutate(ttt2, 
                      total = non_adapted + adapted,
                      perc = round(adapted / (adapted + non_adapted), 2))

clone_counts_2N <- ggplot() + 
  geom_bar(data = ttt, aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', position = position_dodge(0.5), alpha = 0.5) +
  theme_bw() + scale_fill_manual(values = c('gray40','dodgerblue'), name = 'adapted') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("n BCs w/fitness") + 
  xlab("source env") +
  annotate(geom = "text", label = paste0(round(ttt2$perc,2)), x = ttt2$source, y = -20, size = 2) +
  scale_y_continuous(limits = c(-20,600), breaks = seq(0,600,100))

# now tally up with mutation data:
adapteds_wgs <- adapteds[adapteds$wgs == TRUE, ]

mmm <- data.frame(table(adapteds_wgs$is_adapted,
                        adapteds_wgs$Subpool.Environment)) %>% 
  dplyr::arrange(desc(Var1),desc(Freq))

# scrub levels that don't make sense; re-order:
mmm <- mmm[!mmm$Var2 %in% c('none','not_read'),]
mmm$Var2 <- factor(mmm$Var2, levels = paste0(unique(mmm$Var2)))

# supply text annotations:
mmm2 <- reshape2::dcast(mmm, Var2 ~ Var1, value.var = 'Freq')
names(mmm2) <- c('source','non_adapted','adapted')
mmm2 <- dplyr::mutate(mmm2,
                      total = non_adapted + adapted,
                      perc = round(adapted / (adapted + non_adapted),2))

wgs_counts_2N <- ggplot() + 
  geom_bar(data = mmm, aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', position = position_dodge(0.5), alpha = 0.5) +
  theme_bw() + scale_fill_manual(values = c('gray40','dodgerblue'), name = 'adapted') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("n BCs w/wgs data") + 
  xlab("source env") +
  annotate(geom = "text", label = paste0(round(mmm2$perc,2)), x = mmm2$source, y = -5, size = 2) +
  scale_y_continuous(limits = c(-5,120), breaks = seq(0,120,20))

# generate plottable object:
clone_counts_bar1 <- ggarrange(plotlist = list(clone_counts_2N, wgs_counts_2N), labels = c('a','b'), 
                               hjust = 0,
                               vjust = 1,
                               align = 'hv',
                               widths = c(1,0.66), common.legend = T)


# generate fitness clusters:
x <- adapteds
the_cols <- focal_cols[!focal_cols %in% c('Diverse.BC',
                                          'Environment.BC',
                                          'maha_dist',
                                          'dist_pval',
                                          'is_adapted',
                                          'neutral_set',
                                          'Which.Subpools',
                                          grep('_err$', focal_cols, value = T),
                                          grep(arguments$exclude, focal_cols, value = T))]

plot_fitness <- function(x, the_cols) {
  
  x <- 
    x %>%
    dplyr::filter(is_adapted == TRUE) %>%
    dplyr::select(the_cols, wgs) %>%
    dplyr::mutate(CLM.iva_s = CLM.iva_s / 10,
                  FLC4.iva_s = FLC4.iva_s / 10) %>%
    tidyr::gather(key = 'bfa_env', value = 'fitness', -Full.BC, -Subpool.Environment, -wgs) %>%
    dplyr::filter(!Subpool.Environment %in% c('none','not read'))
  
  fit_plot_1 <- ggplot(x, aes(x = bfa_env, y = fitness / as.double(arguments$gens), group = Full.BC, col = Subpool.Environment)) +
    facet_grid(Subpool.Environment ~ wgs) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 6),
          legend.position = 'none') +
    geom_line(alpha = 0.33) +
    ylab('s per gen') +
    geom_hline(yintercept = 0, col = 'black') +
    scale_y_continuous(limits = c(-0.2,0.1), breaks = seq(-0.2,0.1,0.025))
  
}