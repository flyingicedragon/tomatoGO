#' Submit revigo job
#'
#' @param go_result GO result with tab seperation.
#' @param revigo_out Output path.
#'
#' @export
revigo <- function(go_result, revigo_out) {
  revigo_cmd <- paste(
    "src/revigo.py",
    "--go", go_result,
    "--data",
    "--out", revigo_out
  )
  system(revigo_cmd)
}

#' @title Read revigo table
#'
#' @param path path to revigo table
#'
#' @importFrom dplyr filter mutate
#' @importFrom readr read_tsv
#' @importFrom rlang .data
#' @return revigo tibble
#' @export
read_revigo <- function(path) {
  result <- read_tsv(path) %>%
    filter(.data[["eliminated"]] == 0) %>%
    mutate(
      plot_X = as.numeric(.data[["plot_X"]]),
      plot_Y = as.numeric(.data[["plot_Y"]])
    )
  result
}

#' @title Plot revigo results
#'
#' @param result filtered result of revigo
#' @inheritDotParams plot_revigo_scatter_one -result_one
#' @importFrom dplyr filter
#' @return ggplot2 object list
#' @export
plot_revigo_scatter <- function(result, ...) {
  pic <- list()
  for (item in c(
    "Biological_Process", "Cellular_Component", "Molecular_Function"
  )) {
    result_one <- filter(result, `%namespace` == item)
    pic[[item]] <- plot_revigo_scatter_one(result_one, ...)
  }
  pic
}

#' @title Plot revigo results per item
#'
#' @param result_one one result item
#' @param top top dispensability number for show
#' @param ... label show condition
#'
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_gradientn
#' @importFrom ggplot2 scale_size theme_bw theme labs element_blank
#' @importFrom ggplot2 element_text xlim ylim
#' @importFrom dplyr top_n
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
