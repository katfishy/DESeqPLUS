library(testthat)
library(DESeqPLUS)
library(DESeq2)
library(airway)
data(airway)
dds <- DESeqDataSet(airway, design = ~ dex)

test_that("qcPlot returns a list of ggplot objects", {
  plots <- qcPlot(dds, min_counts = 10)
  
  # Test outputs
  expect_type(plots, "list")
  expect_length(plots, 2)
  expect_s3_class(plots$library_plot, "gg")
  expect_s3_class(plots$box_plot, "gg")
  
  # Test error handling
  expect_error(qcPlot(NULL), "Function input must be in DESEqDataSet format.")
})