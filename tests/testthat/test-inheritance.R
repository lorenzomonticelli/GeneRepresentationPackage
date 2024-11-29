# test-inheritance

test_that("Subclasses inherit from Gene", {
  pc_structure <- GRanges(seqnames = "chr17", ranges = IRanges(start = 7668402, end = 7687550), strand = "+")
  pc_gene <- createProteinCodingGene(
    "ENSG000001", "TP53", "Tumor Protein", "Tumor suppressor gene",
    pc_structure, "P04637", "MVLSPADKTNVKAAWG"
  )
  expect_true(is(pc_gene, "Gene"))

  lnc_structure <- GRanges(seqnames = "chr11", ranges = IRanges(start = 65497758, end = 65503100), strand = "+")
  lnc_gene <- createLncRNAGene(
    "ENSG000002", "MALAT1", "Metastasis Associated Lung Adenocarcinoma Transcript 1",
    "lncRNA associated with metastasis", lnc_structure, "NONHSAG000002", "AUGCUACGUGA"
  )
  expect_true(is(lnc_gene, "Gene"))

  mir_structure <- GRanges(seqnames = "chr21", ranges = IRanges(start = 25575833, end = 25575902), strand = "-")
  mir_gene <- createMiRNAGene(
    "ENSG000003", "MIR21", "microRNA 21", "A microRNA associated with cancer",
    mir_structure, "MI0000077", "AGCUUA"
  )
  expect_true(is(mir_gene, "Gene"))
})
