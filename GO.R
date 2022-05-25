library(readr)
library(dplyr)
library(this.path)
library(clusterProfiler)

# import GO annotation -----

go_annotation_all <- read_csv(
  here("tomato_go.csv")
)
go_annotation <- split(go_annotation_all, with(go_annotation, level))
go_info <- read_tsv(
  here("go-basic.tb")
)
type_vector <- c(
  "biological_process",
  "cellular_component",
  "molecular_function"
)

# GO

#' @title split GO enricher
#' @param genes vector to do GO analysis
#' @param GO types
#' @return splitted GO enricher
go_split_enricher <- function(gene_list, type) {
  go_result_split <- enricher(
    gene_list,
    pAdjustMethod = "none",
    TERM2GENE = go_annotation[[type]][c(2, 1)],
    TERM2NAME = go_info[1:2]
  )
  return(go_result_split)
}

#' @title GO analysis per type
#' @param gene_list genes vector
#' @param sample_name sample name used in result file
go_enricher <- function(gene_list, sample_name) {
  go_result <- list()
  for (ii in type_vector) {
    go_result[[ii]] <- go_split_enricher(gene_list, ii)
  }
  barplot_save <- function(go_result, go_type) {
    img <- barplot(go_result, showCategory = 30)
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
    barplot_save(go_result[[ii]], ii)
  }
  return(go_result)
}

#' @title GO analysis for all
#' @param gene_list genes vector
#' @param sample_name sample name used in result file
go_enricher_all <- function(gene_list, sample_name) {
  go_result <- go_split_enricher(gene_list, sample_name)
  ## bar plot
  barplot_save <- function(go_result) {
    img <- barplot(go_result, showCategory = 30)
    ggsave(
      paste0(sample_name, ".pdf"),
      plot = img,
      width = 12,
      height = 6
    )
  }
  barplot_save(go_result)
  return(go_result)
}

#' @title get GO genes
#' @param GO result
#' @return GO result data
get_genes <- function(go_result) {
  type_vector <- c(
    "biological_process",
    "cellular_component",
    "molecular_function"
  )

  get_go_list <- function(go_type_result) {
    go_type_result <- tibble::as_tibble(go_type_result)[, c(1, 2, 8)]
    return(go_type_result)
  }

  go_result_list <- list()
  for (ii in type_vector) {
    go_result_list[[ii]] <- get_go_list(go_result[[ii]])
    go_result_list[[ii]] <- go_result_list[[ii]] %>%
      separate_rows(geneID, sep = "/")
    go_result_list[[ii]] <- dplyr::mutate(go_result_list[[ii]], type = ii)
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

#' @title get GO genes
#' @param gene_list genes list for GO analysis
#' @param example_name name to used in file
#' @return tibble containing GO genes
get_go_genes <- function(gene_list, example_name) {
  gene_go <- go_enricher(gene_list, example_name) %>%
    get_genes()
  return(gene_go)
}