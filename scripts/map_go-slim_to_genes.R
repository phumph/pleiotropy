#!/usr/bin/env Rscript

# map_go-slim_to_gene.R

# Sets up and maps GO-slim categories per gene in mutation list.


suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))


MUTATIONS_FILE <- "data/mutation_data/mutations_by_bc.csv"
GO_OUTFILE <- "data/mutation_data/go-slim_results.txt"
GO_PARSED_OUTFILE <- "data/mutation_data/go-slim_by_gene_parsed.txt"


generate_mut_output <- function(mut_df, outfile) {
  
  mut_df %>%
    dplyr::select(GENE) %>%
    unique() %>%
    dplyr::mutate(QUERY_NAME = GENE) ->
    muts_for_go
  
  # manually updated query name to enable retrieval of data for all genes
  # from yeastgenome.org
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "HAP1"] <- "YLR256W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "MIF2"] <- "YKL089W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "MLP1"] <- "YKR095W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "SLM2"] <- "YNL047C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "SNM1"] <- "YDR478W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "TAF1"] <- "YGR274C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "TEX1"] <- "YNL253W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "TIP1"] <- "YBR067C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "CIS1"] <- "YLR346C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "ANT1"] <- "YPR128C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "CTH1"] <- "YDR151C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "HAP3"] <- "YBL021C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "POT1"] <- "YIL160C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "RET1"] <- "YOR207C"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "SDL1"] <- "YIL167W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "SHE1"] <- "YBL031W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "TUF1"] <- "YOR187W"
  muts_for_go$QUERY_NAME[muts_for_go$GENE == "YOL013W-B"] <- "S000007252"
  
  muts_for_go$QUERY_NAME %>% 
    sort() ->
    gene_list
  
  gene_list %>%
    readr::write_lines(path = outfile)
  
  return(gene_list)
}


split_row <- function(go, genes) {
  genes %>% 
    strsplit(",") %>% unlist() -> 
    genes
  sapply(genes, function(x) gsub(" ", "", x)) ->
    genes
  df <- data.frame(GOID = go, GENE = genes, stringsAsFactors = FALSE)
  row.names(df) <- NULL
  return(df)
}


parse_go_data <- function(go_df, GOID = "GOID", GENEID = "ANNOTATED_GENES") {
  
  go_df %>%
    dplyr::select(GOID, ANNOTATED_GENES) ->
    go_slim_genes
  
  lapply(split(go_slim_genes,
               go_slim_genes[[GOID]]), function(x) {
                 split_row(go = x[[GOID]],
                           genes = x[[GENEID]])
               }
  ) %>% do.call(rbind, .) ->
    go_slim_long
  
  row.names(go_slim_long) <- NULL
  
  go_slim_long %>%
    dplyr::mutate(hit = 1) %>%
    tidyr::pivot_wider(names_from = GOID, values_from = hit, values_fill = list("hit" = 0)) ->
    go_slim_wide
  
  return(list("go_long" = go_slim_long,
              "go_wide" = go_slim_wide))
}


main <- function() {
  
  mutation_data <- read.table(MUTATIONS_FILE,
                              sep = ",",
                              header = TRUE,
                              stringsAsFactors = FALSE)
  
  gene_list <- generate_mut_output(mutation_data, outfile = "data/mutation_data/gene_list.txt")
  
  go_slim_dat <- read.table(GO_OUTFILE,
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE)
  
  go_dat_parsed <- parse_go_data(go_slim_dat)
  
  go_dat_parsed$go_long %>%
    readr::write_csv(path = GO_PARSED_OUTFILE)
}

main()
