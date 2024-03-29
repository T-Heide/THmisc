#' GGplot layer to plot a violing plot
#'
#' @param color Color of violin plot (default: light gray).
#' @param fill Color of overlayed boxplot (default: dark gray).
#' @param ... unused
#'
#' @return A list containg ggplot layers to create a violin plot with a overlayed boxplot and mean and median points.
#' @export
#'
pretty_violin =
  function(color="gray85", fill="gray10", ...) {
    list(
      ggplot2::geom_violin(fill=color, scale = "width", width=0.8),
      ggplot2::geom_boxplot(fill=fill, color=fill, outlier.shape = NA, width=0.2),
      ggplot2::stat_summary(fun = "mean", aes(shape="Mean"), color=color),
      ggplot2::stat_summary(fun = "median", aes(shape="Median"), color=color, size=2),
      ggplot2::labs(shape=""),
      ggplot2::scale_shape_manual(values=c(20, 95), breaks=c("Mean","Median")),
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color="gray10", fill=NA, linetype=0))),
      cowplot::theme_cowplot()
    )
  }
