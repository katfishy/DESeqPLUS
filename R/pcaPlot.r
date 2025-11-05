#' Principal Component Analysis (PCA) on Normalized RNA-seq Data
#'
#' This function produces a PCA plot to explore sample relationships and
#' detect possible outliers of the DESeqDataSet file.
#'
#' @param dds : DESeqDataSet
#' @param ntop : Number of most variable genes to use for PCA (default is 500)
#' @param title : Optional title for plot
#'
#' @return a ggplot2 PCA plot
#'
#' @examples library(DESeq2)
#' library(airway)
#' data(airway)
#' dds <- DESeqDataSet(airway, design = ~ dex)
#' pcaPlot(dds)
#'
#' @references
#' Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and
#' dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014).
#' https://doi.org/10.1186/s13059-014-0550-8
#'
#' Morgan M, Obenchain V, Hester J, Pagès H (2025). SummarizedExperiment:
#' A container (S4 class) for matrix-like assays.
#' doi:10.18129/B9.bioc.SummarizedExperiment, R package version 1.40.0,
#' https://bioconductor.org/packages/SummarizedExperiment.
#'
#' OpenAI. (2025). ChatGPT (GPT-5) [Large language model] https://chat.openai.com/
#'
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot aes_string geom_point theme_minimal labs theme element_text
#' @importFrom stats prcomp
#' @export
pcaPlot <- function(dds, ntop=500, title="PCA Plot") {
  if (inherits(dds, "DESeqDataSet")) {
    # Normalize data
    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
    mat <- SummarizedExperiment::assay(vsd)
    metadata <- as.data.frame(SummarizedExperiment::colData(vsd))
  } else {
    stop("Function input must be in DESEqDataSet format.")
  }

  vars <- apply(mat, 1, var, na.rm=TRUE)
  top_genes <- order(vars, decreasing=TRUE)[seq_len(min(ntop, length(vars)))]
  mat_top <- t(mat[top_genes, ])

  # PCA
  pca <- stats::prcomp(mat_top, scale.=TRUE)
  pca_data <- data.frame(pca$x, metadata)

  # Percent Variance
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

  color_col <- names(metadata)[1]
  if (is.null(color_col)) {
    color_col <- "Sample"
    metadata$Sample <- colnames(mat)
    pca_data <- data.frame(pca$x, metadata)
  }

  # Plot
  plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1,
                                                 y = PC2,
                                                 color = pca_data[[color_col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = color_col
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(plot)
}
