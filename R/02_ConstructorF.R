#' Validate GRanges Structure
#'
#' Checks whether the provided structure is a valid GRanges object.
#'
#' @importFrom GenomicRanges GRanges
#' @param structure A GRanges object containing gene structure information.
#' @return Throws an error if the structure is not a GRanges object.
#' @keywords internal
validateGRanges <- function(structure) {
    if (!inherits(structure, "GRanges")) {
        stop("The 'structure' slot must be a GRanges object.")
    }
}

#' Create a Protein Coding Gene Object
#'
#' Constructor function for creating a ProteinCodingGene object.
#'
#' @param id Character. Gene identifier.
#' @param symbol Character. Gene symbol.
#' @param name Character. Gene name.
#' @param description Character. Gene description.
#' @param structure GRanges. Gene structure data.
#' @param protein_id Character. Protein identifier.
#' @param protein_sequence Character. Protein sequence.
#' @return An instance of the ProteinCodingGene class.
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
#' @export
createProteinCodingGene <- function(id, symbol, name, description, 
    structure, protein_id, protein_sequence) {
    validateGRanges(structure)
    new("ProteinCodingGene",
        id = id, symbol = symbol, name = name, description = description,
        structure = structure, protein_id = protein_id, 
        protein_sequence = protein_sequence
    )
}

#' Create a Long Non-Coding RNA Gene Object
#'
#' Constructor function for creating an lncRNAGene object.
#'
#' @param id Character. Gene identifier.
#' @param symbol Character. Gene symbol.
#' @param name Character. Gene name.
#' @param description Character. Gene description.
#' @param structure GRanges. Gene structure data.
#' @param lncRNA_id Character. lncRNA identifier.
#' @param RNA_sequence Character. RNA sequence.
#' @return An instance of the lncRNAGene class.
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
#' @export
createLncRNAGene <- function(id, symbol, name, description, 
    structure, lncRNA_id, RNA_sequence) {
    validateGRanges(structure)
    new("lncRNAGene",
        id = id, symbol = symbol, name = name, description = description,
        structure = structure, lncRNA_id = lncRNA_id, 
        RNA_sequence = RNA_sequence
    )
}

#' Create a microRNA Gene Object
#'
#' Constructor function for creating a miRNAGene object.
#'
#' @param id Character. Gene identifier.
#' @param symbol Character. Gene symbol.
#' @param name Character. Gene name.
#' @param description Character. Gene description.
#' @param structure GRanges. Gene structure data.
#' @param miRNA_id Character. microRNA identifier.
#' @param seed_sequence Character. Seed sequence of the microRNA.
#' @return An instance of the miRNAGene class.
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
#' @export
createMiRNAGene <- function(id, symbol, name, description, 
    structure, miRNA_id, seed_sequence) {
    validateGRanges(structure)
    new("miRNAGene",
        id = id, symbol = symbol, name = name, description = description,
        structure = structure, miRNA_id = miRNA_id, 
        seed_sequence = seed_sequence
    )
}
