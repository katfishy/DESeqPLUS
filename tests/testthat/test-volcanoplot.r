library(testthat)
library(DESeqPLUS)
library(DESeq2)
library(airway)

data(airway)
dds <- DESeqDataSet(airway, design = ~ dex)
dds <- DESeq(dds)
res <- results(dds)

test_that("volcanoPlot works correctly", {
  plot <- volcanoPlot(res, lfc_threshold = 1, pval_threshold = 0.1)
  
  # Test output
  expect_s3_class(plot, "gg")
  
  # Test error handling
  expect_error(volcanoPlot(NULL), "Function input must be in DESeqResults format.")
})