library(tidyverse)
library(clusterProfiler)

go_enricher <- function(gene_list, sample_name) {
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
}
