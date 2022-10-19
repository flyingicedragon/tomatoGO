#' Do GO analysis for tomato with all type
#' @param go_type GO types
#' @param ... other parameters for enricher
#' @inheritDotParams clusterProfiler::enricher gene pvalueCutoff:qvalueCutoff
#' @importFrom clusterProfiler enricher
#' @return GO retults
#' @export
tomato_go <- function(gene, go_type = "all", ...) {
  if (go_type == "all") {
    go_gson <- tomato_all
  } else if (go_type == "bp") {
    go_gson <- tomato_bp
  } else if (go_type == "cc") {
    go_gson <- tomato_cc
  } else if (go_type == "mf") {
    go_gson <- tomato_mf
  }
  go_result <- enricher(
    gene,
    gson = go_gson,
    ...
  )
  go_result
}

#' @title Get GO genes
#' @param go_result GO result
#' @param all GO result contains all GO type or not
#' @importFrom tidyr separate_rows
#' @importFrom dplyr "%>%"
#' @return GO result data tibble
#' @export
get_genes <- function(go_result, all = TRUE) {
  if (all) {
    go_result <- tibble::as_tibble(go_result) %>%
      separate_rows(geneID, sep = "/") %>%
      as.matrix()
  } else {
    get_go_list <- function(go_type_result) {
      go_type_result <- tibble::as_tibble(go_type_result)
      go_type_result
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
  }
  go_result <- tibble::as_tibble(go_result)
  go_result
}

#' @title Get GO results
#' @param result GO result
#' @param value value names for export. eg. pvalue, p.adjust, qvalue, Count
#' @importFrom dplyr select
#' @return tibble containing results
get_go_results <- function(result, value = "pvalue") {
  go_result <- tibble::as_tibble(result)
  go_result <- select(go_result, "ID", !!value)
  go_result
}

#' @title Export GO result with p-value
#'
#' @param result_all GO result
#' @param path path to save result
#' @inheritDotParams get_go_results value
#'
#' @importFrom readr write_tsv
#' @export
export_go_results <- function(result_all, path, ...) {
  result <- get_go_results(result_all, ...)
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