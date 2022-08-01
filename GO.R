library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(this.path)
library(clusterProfiler)

# import GO annotation -----

go_annotation_all <- read_tsv(
  here("tomato_go.tsv"),
  col_names = c("id", "go", "level")
)
go_annotation <- split(go_annotation_all, with(go_annotation_all, level))
go_info <- read_tsv(
  here("go-basic.tb")
)
type_vector <- c(
  "biological_process",
  "cellular_component",
  "molecular_function"
)

# GO ----

#' @title split GO enricher
#' @param genes vector to do GO analysis
#' @param GO types
#' @return splitted GO enricher
go_split_enricher <- function(gene_list, type = FALSE) {
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
  go_result <- enricher(
    gene_list,
    pAdjustMethod = "none",
    TERM2GENE = go_annotation_all[c(2, 1)],
    TERM2NAME = go_info[1:2]
  )
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

#' @title get GO results
#' @param genes vector contains genes
#' @param value value names for export. eg. pvalue, p.adjust, qvalue, Count
#' @return tibble containing results
get_go_results <- function(genes, value = "pvalue") {
  go_result <- enricher(
    genes,
    pAdjustMethod = "none",
    TERM2GENE = go_annotation_all[c(2, 1)],
    TERM2NAME = go_info[1:2]
  ) %>%
    tibble::as_tibble() %>%
    select("ID", !!value)
  go_result
}

# revigo -----

#' @title Export GO result with p-value
#'
#' @param genes genes list
export_go_results <- function(genes, path) {
  result <- get_go_results(genes[["name"]])
  write_tsv(
    result,
    path,
    col.names = FALSE
  )
  shell_str <- paste0(
    "java -jar /mnt/hdd0/icedragon/packages/",
    "RevigoStandalone_2015-02-17_beta/RevigoStandalone.jar ",
    path, " --taxId=4081"
  )
  print(shell_str)
}

#' @title Read revigo table
#'
#' @return revigo tibble
read_revigo <- function(path) {
  result <- read_tsv(path) %>%
    dplyr::filter(eliminated == 0) %>%
    mutate(plot_X = as.numeric(plot_X), plot_Y = as.numeric(plot_Y))
  result
}

#' @title Plot revigo results
#'
#' @param result filtered result of revigo
plot_revigo_scatter <- function(result) {
  pic <- list()
  for (item in c(
    "Biological_Process", "Cellular_Component", "Molecular_Function"
  )) {
    result_one <- dplyr::filter(result, `%namespace` == item)
    pic[[item]] <- plot_revigo_scatter_one(result_one)
  }
  pic
}

#' @title Plot revigo results per item
#'
#' @param result_one one result item
#'
#' @return ggplot2 object
plot_revigo_scatter_one <- function(result_one) {
  p1 <- ggplot(data = result_one)
  p1 <- p1 + geom_point(aes(
    plot_X, plot_Y,
    size = plot_size, color = `log10_p-value`,
    alpha = I(0.6)
  ))
  p1 <- p1 + scale_colour_gradientn(
    colors = c("red", "yellow", "green"),
    limits = c(min(result_one[["log10_p-value"]]), 0)
  )
  p1 <- p1 + scale_size(range = c(5, 30)) + theme_bw()
  ex <- dplyr::filter(
    result_one,
    dispensability < 0.1 |
      label == 1 # |
    # `log10_p-value` < -6
  )
  p1 <- p1 + geom_text(
    data = ex,
    aes(plot_X, plot_Y, label = description),
    color = I(alpha("black", 0.85)),
    size = 5
  )
  p1 <- p1 + labs(y = "semantic space x", x = "semantic space y")
  p1 <- p1 + theme(legend.key = element_blank())
  p1 <- p1 + theme(text = element_text(size = 18))
  result_x <- max(result_one$plot_X) - min(result_one$plot_X)
  result_y <- max(result_one$plot_Y) - min(result_one$plot_Y)
  p1 <- p1 + xlim(
    min(result_one$plot_X) - result_x / 10,
    max(result_one$plot_X) + result_x / 10
  )
  p1 <- p1 + ylim(
    min(result_one$plot_Y) - result_y / 10,
    max(result_one$plot_Y) + result_y / 10
  )
  p1
}