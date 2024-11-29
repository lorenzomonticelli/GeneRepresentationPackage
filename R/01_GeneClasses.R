#' Gene Class
#'
#' A virtual S4 class to represent the generalities of genes.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import methods
#' @slot id Character. The gene identifier (e.g., Ensembl or NCBI gene ID).
#' @slot symbol Character. The HUGO symbol of the gene.
#' @slot name Character. The full name of the gene.
#' @slot description Character. The description of the gene.
#' @slot structure GRanges. The genomic structure of the gene 
#' (chromosome, start, end, strand).
#' @param object An object of class Gene.
#' @param value A value to be assigned to the slot.
#' @examples
#' library(GenomicRanges)
#' # Define a generic genomic structure
#' structure <- GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges(start = 1000, end = 5000),
#'   strand = "+"
#' )
#' # Create a generic gene object
#' gene <- new(
#'   "Gene",
#'   id = "GENE000001",
#'   symbol = "GENE1",
#'   name = "Generic Gene",
#'   description = "This is a generic example of a gene.",
#'   structure = structure
#' )
#' print(gene)
#' @keywords classes
#' @export
setClass("Gene",
    slots = list(
        id = "character",
        symbol = "character",
        name = "character",
        description = "character",
        structure = "GRanges"
    )
)

#' Protein Coding Gene Class
#'
#' An S4 class representing protein-coding genes,
#' inheriting general properties from the Gene class.
#'
#' @slot protein_id Character. Protein identifier.
#' @slot protein_sequence Character. Protein sequence.
#' @examples
#' library(GenomicRanges)
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
#' print(tp53)
#' @keywords classes
#' @export
setClass("ProteinCodingGene",
    contains = "Gene",
    slots = list(
        protein_id = "character",
        protein_sequence = "character"
    )
)

#' Long Non-Coding RNA (lncRNA) Gene Class
#'
#' An S4 class representing long non-coding RNA (lncRNA) genes,
#' inheriting properties from the Gene class.
#'
#' @slot lncRNA_id Character. lncRNA identifier.
#' @slot RNA_sequence Character. RNA sequence.
#' @examples
#' library(GenomicRanges)
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
#' print(malat1)
#' @keywords classes
#' @export
setClass("lncRNAGene",
    contains = "Gene",
    slots = list(
        lncRNA_id = "character",
        RNA_sequence = "character"
    )
)

#' microRNA Gene Class
#'
#' An S4 class representing microRNA (miRNA) genes, 
#' inheriting properties from the Gene class.
#'
#' @slot miRNA_id Character. microRNA identifier.
#' @slot seed_sequence Character. Seed sequence of the microRNA.
#' @examples
#' library(GenomicRanges)
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
#' print(mir21)
#' @keywords classes
#' @export
setClass("miRNAGene",
    contains = "Gene",
    slots = list(
        miRNA_id = "character",
        seed_sequence = "character"
    )
)
