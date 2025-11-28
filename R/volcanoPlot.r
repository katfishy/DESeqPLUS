#' Volcano Plot from DESeq2 results
#'
#' This function produces a volcano plot showing log2 fold changes vs -log10 p-values.
#' Points with p-value less than threshold are colored in red, and others are grey.
#'
#' @param res : DESeqResults
#' @param lfc_threshold : Log2 fold change threshold for highlighting (default is 0)
#' @param pval_threshold : Adjusted p-value threshold (default is 0.05)
#' @param title : Optional title for the plot
#'
#' @return a ggplot2 volcano plot
#'
#' @examples library(DESeq2)
#' library(airway)
#' data(airway)
#' dds <- DESeqDataSet(airway, design = ~ dex)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' volcanoPlot(res)
#'
#' @references
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs theme element_text
#' @export
volcanoPlot <- function(res, lfc_threshold=0, pval_threshold=0.05, title="Volcano Plot") {
  if (!("DESeqResults" %in% class(res))) {
    stop("Function input must be in DESeqResults format.")
  }

  df <- data.frame(
    log2FoldChange = res$log2FoldChange,
    pval = res$padj
  )
  df <- df[!is.na(df$pval), ]

  df$color <- ifelse(df$pval < pval_threshold & abs(df$log2FoldChange) > lfc_threshold, "red", "gray")

  # Plot
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = -log10(pval), color = color)) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significant"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  return(plot)
}
