#### plt_heatmap1.R
## Produces heatmap from parsed BFA file
## Last updated: 11 April 2018 by PTH


#### Header ####
library("ggplot2")
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


#### function definitions ####
Inf_to_NA <- function(x){
  # find rows with Inf:
  inf.rows <- apply(x[,-1],1,function(x) sum(is.infinite(x)))
  
  # find submatrix with Infs:
  sm <- x[inf.rows>0 , ]
  
  # traverse rows
  for (r in 1:length(sm[,1])){
    sm[r,-1][is.infinite(as.numeric(as.vector(sm[r,-1])))] <- NA
  }
  
  # substitute Inf-->NA rows back into input df
  x[inf.rows>0,] <- sm
  
  # output modified input df
  return(x)
}

# define custom dist and clust functions for use with heatmap.3
#mydist=function(c) {daisy(Inf_to_NA(c))}
mydist=function(c) {dist(c, method = 'euclidean')}
myclust=function(c) {hclust(c, method='average')}

# function to transform fitness data into quantile (ranks) per bfa.env
sqr_by_env <- function(x){
  # define function to calculate signed quantile rescaling
  sqr <- function(x){
    x2 <- x>0 #take positive
    x3 <- x<0 #take negatives
    # re-scale positives
    x[x2] <- rank(x[x2]) / length(x[x2])
    x[x3] <- (-1) * rank(abs(x[x3])) / length(x[x3])
    return(x)
  }
  
  # apply function over all input columns separately
  x2 <- do.call(cbind, lapply(x, sqr))
  
  # return transformed 
  return(x2)
}

#### Prune input data to focal envs and sources ####
# import data
BFA1 <- read.csv(here("data/dBFA2_hBFA1_with_adapt_11APR2018.csv"))

# remove source==48Hr BCs:
BFA2 <- BFA1[-c(grep('48Hr', BFA1$Subpool.Environment),
                grep('48Hr', BFA1$bfa.env),
                grep('X37C_Stan', BFA1$bfa.env)),]

# split data into two: one with meta-data, the other with fitnesses, then re-join:
bfa_meta <- unique(BFA2[,names(BFA2) %in% c('Full.BC','adapted','neutral','source','ploidy','subpool','autodip')])
BFA3 <- unique(BFA2[,names(BFA2) %in% c('Full.BC','bfa.env','iva.s')])
BFA4 <- reshape2::dcast(BFA3, Full.BC ~ bfa.env, value.var = 'iva.s')

# add back meta-data:
BFA5 <- merge(BFA4,bfa_meta, by = 'Full.BC', sort = F)
row.names(BFA5) <- BFA5$Full.BC

# define bfa.env names:
bfa.envs <- names(BFA5)[2:11]

# transform any 'Inf' values to 'NA'
BFA6 <- Inf_to_NA(BFA5)

# get rid df of rows with too many NAs
NAs <- apply(BFA6[,-1], 1, function(x) sum(is.na(x)))
BFA7 <- BFA6[!NAs>2,]

# for now, assign all existing NAs a value of zero.
BFA7[is.na(BFA7)] <- 0
# need to implement distance and clustering algorithm that doesn't break when including NAs.
# Probably the NA values are super low fitness, which in the drug environments is indistinguishable from neutral (i.e. s = 0).

#### Assign master df for all plots ####
BFA <- droplevels(BFA7) # change this if editing above; 'BFA' works with all below:

#### write pruned df ####
write.csv(BFA, file = here("data/dBFA2_hBFA1_with_adapt_11APR2018_pruned.csv"))

#### Define heatmap color properties ####
meta_cols <- c('source','ploidy','subpool','autodip','neutral')
BFA_meta <- BFA[,meta_cols]

BFA_meta$source  <- factor(BFA_meta$source)
BFA_meta$subpool <- factor(BFA_meta$subpool)
BFA_meta$ploidy  <- factor(BFA_meta$ploidy)
BFA_meta$autodip <- factor(BFA_meta$autodip)
BFA_meta$neutral <- factor(BFA_meta$neutral)

# define neutrals colors:
ncolors <- with(BFA_meta,
                data.frame(neutral = levels(neutral),
                           ncolor = c('#f0f0f0','#636363')))

# define auto-diploids colors:
acolors <- with(BFA_meta,
                data.frame(autodip = levels(autodip),
                           acolors = c('#efedf5','#de2d26')))

# define source-ploidy matrix:
colorz <- data.frame(matrix(c("#9ecae1","#3182bd","#a1d99b","#31a354","#fdae6b","#e6550d","#bcbddc","#756bb1","#fc9272","#de2d26","#fa9fb5","#c51b8a","#c994c7","#dd1c77","#8c96c6","#88419d","#fdcc8a","#fc8d59","#cccccc","#525252","#fc8d62","#ff7f00"), ncol = 11))

# define dummy variable matrix for sources:
mm1 <- data.frame(model.matrix(row.names(BFA_meta) ~ 0 + source, data = BFA_meta))
names(colorz) <- names(mm1)

# go through column-wise and assign colors
col_sorter <- function(x, y, z, whitecol = "#f7f7f7"){
  # x is model matrix
  # y is color mapper (row1 = 1N; row2 = 2N)
  # z is character vector of ploidy data
  # for debugging:
  x <- mm1
  y <- colorz
  z <- BFA_meta$ploidy
  c <- 1
  x2 <- x
  # cycle through model matrix columns
  # assign colors for each source according to ploidy:
  for (c in 1:ncol(x)){
    x2[(x[,c]==1) & (z=='1N'), c] <- paste0(y[1,c])
    x2[(x[,c]==1) & (z=='2N'), c] <- paste0(y[2,c])
    # x2[x[z=='1N',c][]==0,c] <- paste0(y[1,c])
    # x2[x[,c][z=='1N']==1,c] <- paste0(y[2,c])
  }
  
  # assign empty color to all remaining zeros:
  x2[x2==0] <- whitecol
  
  # return modified df
  return(x2)
}

mm2 <- col_sorter(x = mm1, y = colorz, z = BFA_meta$ploidy)

# add additional meta-data columns using cbind:
acolors = c('#efedf5','#de2d26')
ncolors = c('#f0f0f0','#636363')

neutrals <- paste0(BFA_meta$neutral)
neutrals[neutrals==0] <- ncolors[1]
neutrals[neutrals==1] <- ncolors[2]

autodips <- paste0(BFA_meta$autodip)
autodips[autodips==0] <- acolors[1]
autodips[autodips==1] <- acolors[2]

# combine all row-wise meta-data for plotting
mm3 <- cbind(mm2, neutrals, autodips)

# make second plot of meta-data,with paneling:
mm3$Full.BC <- row.names(mm3)
mm4 <- reshape2::melt(mm3, measure.vars = names(mm3)[-length(mm3)], id.vars = 'Full.BC', value.name = 'colour', factorsAsStrings = T)

# merge with ploidy info:
temp.df <- data.frame(Full.BC = row.names(BFA_meta), ploidy = BFA_meta$ploidy)

mm4b <- merge(mm4, temp.df, by = 'Full.BC', sort = F)

mm4b$source_ploidy <- paste0(mm4b$variable,'_',mm4b$ploidy)

myColours <- unique(mm4b[,c('source_ploidy','colour')])$colour
names(myColours) <- unique(mm4b[,c('source_ploidy','colour')])$source_ploidy

# need to re-sort according to cluster'd Full.BC, as for hm1:
mm4$Full.BC <- factor(mm4$Full.BC, levels = paste0(labs))

# the color really needs to be barcode and source and ploidy specific.
mm4b$unique.col <- paste0(mm4b$Full.BC,'_',mm4b$variable,'_',mm4b$ploidy)

# re-level as always:
mm4b$Full.BC <- factor(mm4b$Full.BC, levels = paste0(labs))

# re-define colors 
myColours <- mm4b$colour
names(myColours) <- mm4b$unique.col

# re-define order of source environments (to match BFA envs as best as possible)
mm4b$variable <- factor(mm4b$variable, levels = c('sourceGlyEtOH','sourceYPD','sourceAncestor_YPD','sourcepH3_8','sourceSC','source02M_NaCl','source21C','sourcepH7_3','source37C','sourceCLM','sourceFLC4','autodips','neutrals'))

hmeta1 <- ggplot(mm4b, aes(x = variable, y = Full.BC, fill = unique.col)) + geom_tile() + 
  scale_fill_manual(values = myColours) + ylab("") +
  #facet_wrap(~ variable, ncol = 23, scales = "free") +
  theme(axis.text.y = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("meta-data")
#hmeta1


#### Define color scale of heatmap itself ####
# color palette 11-class PuOr
rampcols <- c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b")
rampbreaks = seq(-1,1,0.1)

#### Define input matrix for heatmap.3 ####
map.in <- BFA[,names(BFA) %in% bfa.envs]
map.in <- Inf_to_NA(map.in)
map.in[is.na(map.in)] <- 0

clust1 <- hclust(dist(map.in))
labs <- clust1$labels[clust1$order]

map.in2 <- data.frame(Full.BC = row.names(map.in), sqr_by_env(map.in))
map.in2$Full.BC <- factor(map.in2$Full.BC, levels = paste0(labs))

#### plotting what signed quantiles look like for all environments ####
map.in3 <- reshape2::melt(map.in2)
map.in3$Full.BC <- factor(map.in3$Full.BC, levels = paste0(labs))

# re-order bfa environments
clust2 <- hclust(dist(t(map.in)))
labs2 <- clust2$labels[clust2$order]

map.in3$variable <- factor(map.in3$variable, levels = paste0(labs2))

# colors for gradient
hm_cols <- c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f0f0f0","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b")

# make plot
hm1 <- ggplot(map.in3, aes(x = variable, y = Full.BC, fill = value)) + geom_tile(alpha = 0.5) + 
  scale_fill_gradientn(colors = hm_cols) +
  #scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 0) + ylab("") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("BFA environment") #+
        #ylab("lineage")

#### Adding dendrogram sub-plot ####
#install.packages("ggdendro")
# dd1 <- ggdendrogram(clust1, rotate = )
# dd1

dhc <- as.dendrogram(clust1)
ddata <- dendro_data(dhc, type = "rectangle")
gtree1 <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + theme_dendro() + theme(panel.border = element_blank())

# pdf(file = here("test_hm_combined.pdf"), width = 7, height = 10)
# ggarrange(plotlist = list(gtree1, hmeta1, hm1), align = 'hv', ncol = 3, common.legend = F, widths = c(0.5,0.5,1))
# dev.off()

#### Plot quantile transform correspondence ####
# re-scale data
sqr1 <- sqr_by_env(map.in)
# re-bind to un-scaled for comparison
sqr2 <- cbind(reshape2::melt(sqr1, value.name = 'quantile'),
              reshape2::melt(BFA[, bfa.envs], value.name = 'fitness'))

# plot
pdf(file = here("sqr_fitness_map.pdf"), width = 6, height = 4)
ggplot(test_in3, aes(x = fitness, y = quantile)) + geom_line() + facet_wrap(~ Var2) + 
  geom_vline(xintercept = 0, col = 'gray60', lty = 'dotted') +
  geom_hline(yintercept = 0, col = 'gray60', lty = 'dotted')
dev.off()

