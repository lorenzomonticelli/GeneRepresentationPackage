#' Calculate the Product Length of a Gene
#'
#' Computes the length of the gene product based on the gene type.
#'
#' @importFrom GenomicRanges GRanges
#' @param gene An instance of a Gene or its subclasses 
#' (e.g., ProteinCodingGene, lncRNAGene, miRNAGene).
#' @return An integer representing the length of the gene's product:
#' - Length of the protein sequence for ProteinCodingGene.
#' - Length of the RNA sequence for lncRNAGene.
#' - Length of the seed sequence for miRNAGene.
#' @examples
#' library(GenomicRanges)
#' # Example for ProteinCodingGene
#' structure <- GRanges(
#'   seqnames = "chr17",
#'   ranges = IRanges(start = 7668402, end = 7687550),
#'   strand = "+"
#' )
#' tp53 <- createProteinCodingGene(
#'   id = "ENSG000001",
#'   symbol = "TP53",
#'   name = "Tumor Protein",
#'   description = "Tumor suppressor gene",
#'   structure = structure,
#'   protein_id = "P04637",
#'   protein_sequence = "MVLSPADKTNVKAAWG"
#' )
#' lengthProduct(tp53)
#'
#' # Example for lncRNAGene
#' structure <- GRanges(
#'   seqnames = "chr11",
#'   ranges = IRanges(start = 65497758, end = 65503100),
#'   strand = "+"
#' )
#' malat1 <- createLncRNAGene(
#'   id = "ENSG000002",
#'   symbol = "MALAT1",
#'   name = "Metastasis Associated Lung Adenocarcinoma Transcript 1",
#'   description = "lncRNA associated with metastasis",
#'   structure = structure,
#'   lncRNA_id = "NONHSAG000002",
#'   RNA_sequence = "AUGCUACGUGA"
#' )
#' lengthProduct(malat1)
#'
#' # Example for miRNAGene
#' structure <- GRanges(
#'   seqnames = "chr21",
#'   ranges = IRanges(start = 25575833, end = 25575902),
#'   strand = "-"
#' )
#' mir21 <- createMiRNAGene(
#'   id = "ENSG000003",
#'   symbol = "MIR21",
#'   name = "microRNA 21",
#'   description = "A microRNA associated with cancer",
#'   structure = structure,
#'   miRNA_id = "MI0000077",
#'   seed_sequence = "AGCUUA"
#' )
#' lengthProduct(mir21)
#' @export
setGeneric("lengthProduct", function(gene) {
    standardGeneric("lengthProduct")
})

#' @rdname lengthProduct
#' @aliases lengthProduct,ProteinCodingGene-method
#' @export
setMethod("lengthProduct", "ProteinCodingGene", function(gene) {
    nchar(gene@protein_sequence)
})

#' @rdname lengthProduct
#' @aliases lengthProduct,lncRNAGene-method
#' @export
setMethod("lengthProduct", "lncRNAGene", function(gene) {
    nchar(gene@RNA_sequence)
})

#' @rdname lengthProduct
#' @aliases lengthProduct,miRNAGene-method
#' @export
setMethod("lengthProduct", "miRNAGene", function(gene) {
    nchar(gene@seed_sequence)
})
