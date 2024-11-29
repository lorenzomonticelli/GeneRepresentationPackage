# test-accessor
library(testthat)
library(GenomicRanges)

# Tests for symbol accessors

test_that("symbol accessor works for ProteinCodingGene", {
  structure <- GRanges(seqnames = "chr17", ranges = IRanges(start = 7668402, end = 7687550), strand = "+")
  gene <- createProteinCodingGene(
    "ENSG000001", "TP53", "Tumor Protein", "Tumor suppressor gene",
    structure, "P04637", "MVLSPADKTNVKAAWG"
  )
  expect_equal(symbol(gene), "TP53")

  symbol(gene) <- "TP53_updated"
  expect_equal(symbol(gene), "TP53_updated")
})

test_that("symbol accessor works for lncRNAGene", {
  structure <- GRanges(seqnames = "chr11", ranges = IRanges(start = 65497758, end = 65503100), strand = "+")
  gene <- createLncRNAGene(
    "ENSG000002", "MALAT1", "Metastasis Associated Lung Adenocarcinoma Transcript 1",
    "lncRNA associated with metastasis", structure, "NONHSAG000002", "AUGCUACGUGA"
  )
  expect_equal(symbol(gene), "MALAT1")

  symbol(gene) <- "MALAT1_updated"
  expect_equal(symbol(gene), "MALAT1_updated")
})

test_that("symbol accessor works for miRNAGene", {
  structure <- GRanges(seqnames = "chr21", ranges = IRanges(start = 25575833, end = 25575902), strand = "-")
  gene <- createMiRNAGene(
    "ENSG000003", "MIR21", "microRNA 21", "A microRNA associated with cancer",
    structure, "MI0000077", "AGCUUA"
  )
  expect_equal(symbol(gene), "MIR21")

  symbol(gene) <- "MIR21_updated"
  expect_equal(symbol(gene), "MIR21_updated")
})
