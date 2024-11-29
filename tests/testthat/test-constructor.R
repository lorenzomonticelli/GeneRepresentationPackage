# test-constructor

test_that("createProteinCodingGene initializes correctly", {
  structure <- GRanges(seqnames = "chr17", ranges = IRanges(start = 7668402, end = 7687550), strand = "+")
  gene <- createProteinCodingGene(
    "ENSG000001", "TP53", "Tumor Protein", "Tumor suppressor gene",
    structure, "P04637", "MVLSPADKTNVKAAWG"
  )
  expect_s4_class(gene, "ProteinCodingGene")
  expect_equal(gene@symbol, "TP53")
  expect_equal(gene@protein_sequence, "MVLSPADKTNVKAAWG")
})

test_that("createLncRNAGene initializes correctly", {
  structure <- GRanges(seqnames = "chr11", ranges = IRanges(start = 65497758, end = 65503100), strand = "+")
  gene <- createLncRNAGene(
    "ENSG000002", "MALAT1", "Metastasis Associated Lung Adenocarcinoma Transcript 1",
    "lncRNA associated with metastasis", structure, "NONHSAG000002", "AUGCUACGUGA"
  )
  expect_s4_class(gene, "lncRNAGene")
  expect_equal(gene@symbol, "MALAT1")
  expect_equal(gene@RNA_sequence, "AUGCUACGUGA")
})

test_that("createMiRNAGene initializes correctly", {
  structure <- GRanges(seqnames = "chr21", ranges = IRanges(start = 25575833, end = 25575902), strand = "-")
  gene <- createMiRNAGene(
    "ENSG000003", "MIR21", "microRNA 21", "A microRNA associated with cancer",
    structure, "MI0000077", "AGCUUA"
  )
  expect_s4_class(gene, "miRNAGene")
  expect_equal(gene@symbol, "MIR21")
  expect_equal(gene@seed_sequence, "AGCUUA")
})
