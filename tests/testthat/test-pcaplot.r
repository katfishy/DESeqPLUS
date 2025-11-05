library(testthat)
library(DESeqPLUS)
library(DESeq2)
library(airway)

data(airway)
dds <- DESeqDataSet(airway, design = ~ dex)

test_that("pcaPlot works correctly", {
  plot <- pcaPlot(dds, ntop = 100)
  
  # Test output
  expect_s3_class(plot, "gg")
  
  # Test error handling
  expect_error(pcaPlot(NULL), "Function input must be in DESEqDataSet format.")
})