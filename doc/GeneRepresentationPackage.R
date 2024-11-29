## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GenomicRanges)
library(GeneRepresentationPackage)

## -----------------------------------------------------------------------------
# Define a GRanges object for genomic structure
tp53_structure <- GRanges(
    seqnames = "chr17",
    ranges = IRanges(start = 7668402, end = 7687550),
    strand = "+"
)

# Create a ProteinCodingGene object
tp53 <- createProteinCodingGene(
    id = "ENSG000001",
    symbol = "TP53",
    name = "Tumor Protein",
    description = "Tumor suppressor gene",
    structure = tp53_structure,
    protein_id = "P04637",
    protein_sequence = "MVLSPADKTNVKAAWG"
)
print(tp53)

## -----------------------------------------------------------------------------
# Define a GRanges object for genomic structure
malat1_structure <- GRanges(
    seqnames = "chr11",
    ranges = IRanges(start = 65497758, end = 65503100),
    strand = "+"
)

# Create a lncRNAGene object
malat1 <- createLncRNAGene(
    id = "ENSG000002",
    symbol = "MALAT1",
    name = "Metastasis Associated Lung Adenocarcinoma Transcript 1",
    description = "lncRNA associated with metastasis",
    structure = malat1_structure,
    lncRNA_id = "NONHSAG000002",
    RNA_sequence = "AUGCUACGUGA"
)
print(malat1)

## -----------------------------------------------------------------------------
# Define a GRanges object for genomic structure
mir21_structure <- GRanges(
    seqnames = "chr21",
    ranges = IRanges(start = 25575833, end = 25575902),
    strand = "-"
)

# Create a miRNAGene object
mir21 <- createMiRNAGene(
    id = "ENSG000003",
    symbol = "MIR21",
    name = "microRNA 21",
    description = "A microRNA associated with cancer",
    structure = mir21_structure,
    miRNA_id = "MI0000077",
    seed_sequence = "AGCUUA"
)
print(mir21)

## -----------------------------------------------------------------------------
symbol(tp53)

## -----------------------------------------------------------------------------
symbol(tp53) <- "TP53_updated"
symbol(tp53)

## -----------------------------------------------------------------------------
structure(tp53)

## -----------------------------------------------------------------------------
new_structure <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1000, end = 5000),
    strand = "-"
)
structure(tp53) <- new_structure
structure(tp53)

## -----------------------------------------------------------------------------
lengthProduct(tp53)

## -----------------------------------------------------------------------------
lengthProduct(malat1)

## -----------------------------------------------------------------------------
lengthProduct(mir21)

## ----session-info-------------------------------------------------------------
sessionInfo()

