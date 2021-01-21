library(tidyverse)
library(clusterProfiler)

go_enricher <- function(gene_list, sample_name = "GO") {
  ## import annotation data
  go_annotation <- read_csv("/home/icedragon/DocumentsNotSync/git/tomatoGO/tomato_go.csv")
  go_annotation <- split(go_annotation, with(go_annotation, level))
  go_info <- read_tsv("/home/icedragon/DocumentsNotSync/git/tomatoGO/go-basic.tb")
  type_vector <- c(
    "biological_process",
    "cellular_component",
    "molecular_function"
  )
  ## GO analysis
  go_split_enricher <- function(type) {
    go_result_split <- enricher(
      gene_list,
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
    if (!is.null(go_result) && nrow(go_result) != 0) {
      img <- barplot(go_result, showCategory = 30)
      if (go_type == "molecular_function") {
        img_width <- 12
      } else {
        img_width <- 8
      }
      ggsave(
        paste0(sample_name, "_", go_type, ".pdf"),
        plot = img,
        width = img_width,
        height = 6
      )
    }
  }
  for (ii in type_vector) {
    dotplot_save(go_result[[ii]], ii)
  }
  return(go_result)
}

## get GO genes
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

get_go_genes <- function(gene_list, example_name) {
  gene_go <- go_enricher(gene_list, example_name) %>%
    get_genes()
  return(gene_go)
}
