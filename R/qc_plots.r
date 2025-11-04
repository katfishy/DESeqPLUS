#' Quick Quality Control Plots for RNA-seq data before DESeq2 analysis
#'
#' This function produces a set of QC plots to examine the DESeqDataSet.
#' This includes the total library size per sample and distribution of 
#' counts per sample in log.
#'
#' @param dds : DESeqDataSet
#' @param title : Optional title for plots
#'
#' @return two ggplot2 figures
#'
#' @examples library(DESeq2)
#' library(airway)
#' data(airway)
#' dds <- DESeqDataSet(airway, design = ~ dex)
#' qc_summary(dds)
#' 
#' @importFrom DESeq2 counts
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_minimal labs geom_boxplot theme element_text
#' @export

qc_summary <- function(dds, title="Quality Control Summary") {
  if (inherits(dds, "DESeqDataSet")) {
    counts_data <- counts(dds)
  } else {
    stop("Function input must be in DESEqDataSet format.")
  }
  
  # Total counts per sample
  library_size <- colSums(counts_data)
  df_library <- data.frame(Sample = names(library_size), LibrarySize = library_size)
  
  # Change counts for boxplot
  df_log <- reshape2::melt(log10(counts_data + 1))
  colnames(df_log) <- c("Gene", "Sample", "LogCount")
  
  # Plot 1: Total Library Size
  plot1 <- ggplot2::ggplot(df_library, ggplot2::aes(x=reorder(Sample, LibrarySize), y=LibrarySize)) +
    ggplot2::geom_bar(stat="identity", fill="red") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::labs(title=paste(title, " (Library Size Per Sample)"),
                  x= "Sample",
                  y = "Total Reads")
  
  # Plot 2: Distribution of log10 counts per sample
  plot2 <- ggplot2::ggplot(df_log, ggplot2::aes(x=Sample, y=LogCount, fill = Sample)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::theme_minimal(base_size = 12) +
    # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    #                legend.position = "none") +
    ggplot2::labs(title = paste(title, " (Log10 Counts per Sample)"),
                  x = "Sample",
                  y = "Log10(Counts + 1)")
  
  print(plot1)
  print(plot2)
  
  return(list(library_plot=plot1, box_plot=plot2))
}