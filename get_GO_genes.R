library(tidyverse)
library(clusterProfiler)

go_enricher <- function(gene_list, sample_name) {
  ## import annotation data
  go_annotation <- read_csv("tomato_go.csv")
  go_annotation <- split(go_annotation, with(go_annotation, level))
  go_info <- read_tsv("go-basic.tb")
  type_vector <- c(
    "biological_process",
    "cellular_component",
    "molecular_function"
  )
  ## GO analysis
  go_split_enricher <- function(type) {
    go_result_split <- enricher(
      gene_list[, 1],
      pAdjustMethod = "none",
      TERM2GENE = go_annotation[[type]][c(2, 1)],
      TERM2NAME = go_info[1:2]
    )
    return(go_result_split)
  }
  go_result <- list()
  for (ii in type_vector) {
    go_result[[ii]] <- go_split_enricher(ii)
  }
  ## dot plot
  dotplot_save <- function(go_result, go_type) {
    img <- dotplot(go_result, showCategory = 30)
    if (go_type == "molecular_function") {
      img_width <- 12
    } else {
      img_width <- 8
    }
    ggsave(
      paste0(sample_name, "_", go_type, ".pdf"),
      device = pdf(),
      plot = img,
      width = img_width,
      height = 6
    )
  }
  for (ii in type_vector) {
    dotplot_save(go_result[[ii]], ii)
  }
  return(go_result)
}

read_genes <- function(filename) {
  gene_list <- read.csv(filename)[, 1:2]
  gene_list[, 1] <- str_sub(gene_list[, 1], 1, 14)
  return(gene_list)
}

get_genes <- function(go_result) {
  type_vector <- c(
      "biological_process",
      "cellular_component",
      "molecular_function"
  )
  get_go_list <- function(go_type_result) {
    go_type_result <- as_tibble(go_type_result)[, c(1, 2, 8)]
    return(go_type_result)
  }
  go_result_list <- list()
  for (ii in type_vector) {
    go_result_list[[ii]] <- get_go_list(go_result[[ii]])
    go_result_list[[ii]] <- go_result_list[[ii]] %>%
      separate_rows(geneID, sep = "/")
    go_result_list[[ii]] <- mutate(go_result_list[[ii]], type = ii)
  }
  go_result <- as.matrix(go_result_list[[1]])
  go_result <- rbind(
    go_result,
    go_result_list[[2]] %>% as.matrix()
  ) %>%
    rbind(
      go_result_list[[3]] %>% as.matrix()
    )
  return(go_result)
}

get_go_genes <- function(example_name) {
  gene_list <- read_genes(paste0(example_name, ".csv"))
  gene_go <- go_enricher(gene_list, example_name) %>%
    get_genes()
  gene_go <- merge(
    gene_go,
    gene_list,
    by.x = "geneID",
    by.y = "Accession",
    all.x = T
  )
  names(gene_go) <- c(
    "gene",
    "GO ID",
    "GO term",
    "type",
    "description"
  )
  write.csv(gene_go, paste0(example_name, "_go.csv"), row.names = F)
}
