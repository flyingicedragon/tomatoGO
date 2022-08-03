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

#' @title Split GO enricher
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
  go_result_split
}

#' @title GO analysis per type
#' @param gene_list genes vector
#' @param sample_name sample name used in result file
#' @return GO results list
go_enricher <- function(gene_list, sample_name) {
  go_result <- list()
  for (type in type_vector) {
    go_result[[type]] <- go_split_enricher(gene_list, type)
  }
  go_result
}

#' @title GO analysis for all
#' @param gene_list genes vector
#' @return GO result
go_enricher_all <- function(gene_list) {
  go_result <- enricher(
    gene_list,
    pAdjustMethod = "none",
    TERM2GENE = go_annotation_all[c(2, 1)],
    TERM2NAME = go_info[1:2]
  )
  go_result
}

#' @title Plot GO result
#' @param go_result GO result
#' @param ... other parameters for barplot
#' @return img object
plot_go <- function(go_result, ...) {
  img <- barplot(go_result, showCategory = 30, ...)
  img
}

#' @title Plot GO result for each GO type
#' @param go_result GO results for each GO type
#' @param ... other parameters for barplot
#' @return imgs list
plot_go_each <- function(go_results, ...) {
  imgs <- list()
  for (go_type in type_vector) {
    imgs[go_type] <- plot_go(go_results[go_type], ...)
  }
  imgs
}

#' @title Save plot
#' @param img plot object
#' @param go_type GO type
#' @param path path to save
#' @param sample_name file name to save
plot_save <- function(img, go_type = "all", path = ".", sample_name = "") {
  if (go_type == "molecular_function") {
    img_width <- 12
  } else {
    img_width <- 8
  }
  if (go_type == "all") {
    sample_name <- paste0(sample_name, "_all.pdf")
    sample_name <- file.path(path, sample_name)
    ggsave(
      sample_name,
      device = pdf(),
      plot = img,
      width = img_width,
      height = 6
    )
  } else {
    for (go_type in type_vector) {
      sample_name_one <- paste0(sample_name, "_", go_type, ".pdf")
      sample_name_one <- file.path(path, sample_name_one)
      ggsave(
        sample_name_one,
        device = pdf(),
        plot = img,
        width = img_width,
        height = 6
      )
    }
  }
}

#' @title get GO genes
#' @param go_result GO result
#' @return GO result data matrix
get_genes <- function(go_result, all = FALSE) {
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
  go_result
}

#' @title get GO genes
#' @param gene_list genes list for GO analysis
#' @param example_name name to used in file
#' @return tibble containing GO genes
get_go_genes <- function(gene_list, example_name) {
  gene_go <- go_enricher(gene_list, example_name) %>%
    get_genes()
  gene_go
}

#' @title get GO results
#' @param result GO result all
#' @param value value names for export. eg. pvalue, p.adjust, qvalue, Count
#' @return tibble containing results
get_go_results <- function(result, value = "pvalue") {
  go_result <- tibble::as_tibble(result) %>%
    select("ID", !!value)
  go_result
}

# revigo -----

#' @title Export GO result with p-value
#' @param genes genes vector
#' @param ... other parameters for get_go_results
export_go_results <- function(genes, path, ...) {
  result <- go_enricher_all(genes) %>%
    get_go_results(genes, ...)
  write_tsv(
    result,
    path,
    col_names = FALSE
  )
  shell_str <- paste0(
    "java -jar /mnt/hdd0/icedragon/packages/",
    "RevigoStandalone_2015-02-17_beta/RevigoStandalone.jar ",
    path, " --taxId=4081"
  )
  print(shell_str)
  shell_str
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
plot_revigo_scatter <- function(result, ...) {
  pic <- list()
  for (item in c(
    "Biological_Process", "Cellular_Component", "Molecular_Function"
  )) {
    result_one <- dplyr::filter(result, `%namespace` == item)
    pic[[item]] <- plot_revigo_scatter_one(result_one, ...)
  }
  pic
}

#' @title Plot revigo results per item
#'
#' @param result_one one result item
#' @param condition label show condition
#'
#' @return ggplot2 object
plot_revigo_scatter_one <- function(result_one, top = FALSE, ...) {
  p1 <- ggplot(data = result_one)
  p1 <- p1 + geom_point(aes(
    plot_X, plot_Y,
    size = plot_size, color = `log10_p-value`,
    alpha = I(0.6)
  ))
  p1 <- p1 + scale_colour_gradientn(
    colors = c("red", "red", "yellow", "green"),
    limits = c(min(result_one[["log10_p-value"]]), 0),
    values = scales::rescale(c(min(result_one[["log10_p-value"]]), -10, -3, 0))
  )
  p1 <- p1 + scale_size(range = c(5, 30)) + theme_bw()
  if (top) {
    ex <- top_n(result_one, top, dispensability)
  } else {
    ex <- dplyr::filter(
      result_one,
      ...
      ## dispensability < 0.1 |
      ##   label == 1 |
      ##   `log10_p-value` < -6
    )
  }
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