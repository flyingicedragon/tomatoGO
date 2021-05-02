library(tidyverse)
library(clusterProfiler)

go_enricher <- function(gene_list, sample_name) {
  ## import annotation data
  go_annotation <- read_csv(
    "/home/icedragon/DocumentsNotSync/git/tomatoGO/tomato_go.csv"
  )
  # go_annotation <- split(go_annotation, with(go_annotation, level))
  go_info <- read_tsv(
    "/home/icedragon/DocumentsNotSync/git/tomatoGO/go-basic.tb"
  )
  ## GO analysis
  go_split_enricher <- function() {
    go_result_split <- enricher(
      gene_list[, 1],
      pAdjustMethod = "none",
      TERM2GENE = go_annotation[c(2, 1)],
      TERM2NAME = go_info[1:2]
    )
    return(go_result_split)
  }
  go_result <- go_split_enricher()
  ## dot plot
  dotplot_save <- function(go_result) {
    img <- barplot(go_result, showCategory = 30)
    ggsave(
      paste0(sample_name, ".pdf"),
      plot = img,
      width = 12,
      height = 6
    )
  }
  dotplot_save(go_result)
}
