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
volcanoPlot <- function(res, lfc_threshold=1, pval_threshold=0.05, title="Volcano Plot") {
  if (!("DESeqResults" %in% class(res))) {
    stop("Function input must be in DESeqResults format.")
  }

  df <- data.frame(res)

  df$significance <- "Not Significant"
  df$significance[df$padj < pval_threshold & df$log2FoldChange > lfc_threshold] <- "Up"
  df$significance[df$padj < pval_threshold & df$log2FoldChange < -lfc_threshold] <- "Down"

  # Plot
  plot <- ggplot2::ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted p-value"
    ) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black")

  return(plot)
}
