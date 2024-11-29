#' Gene Accessors
#'
#' Provides accessor methods for retrieving 
#' or updating slots in the Gene class and its subclasses.
#'
#' @import methods
#'
#' @param object An object of class Gene or its subclasses.
#' @param value A value to set for the respective slot.
#'
#' @name GeneAccessors
NULL

#' @rdname GeneAccessors
#' @aliases structure,Gene-method
#' @export
setGeneric("structure", function(object) {
    standardGeneric("structure")
})

#' @rdname GeneAccessors
#' @export
setMethod("structure", "Gene", function(object) {
    object@structure
})

#' @rdname GeneAccessors
#' @export
setGeneric("structure<-", function(object, value) {
    standardGeneric("structure<-")
})

#' @rdname GeneAccessors
#' @export
setReplaceMethod("structure", "Gene", function(object, value) {
    validateGRanges(value)
    object@structure <- value
    validObject(object)
    object
})

#' Get the Symbol of a Gene
#'
#' Retrieves the HUGO symbol of a gene.
#'
#' @inheritParams GeneAccessors
#' @return A character string representing the gene symbol.
#' @examples
#' library(GenomicRanges)
#' gene <- createProteinCodingGene(
#'   id = "ENSG000001",
#'   symbol = "TP53",
#'   name = "Tumor Protein",
#'   description = "Tumor suppressor gene",
#'   structure = GRanges(
#'     seqnames = "chr17",
#'     ranges = IRanges(start = 7668402, end = 7687550),
#'     strand = "+"
#'   ),
#'   protein_id = "P04637",
#'   protein_sequence = "MVLSPADKTNVKAAWG"
#' )
#' symbol(gene)
#' @rdname GeneAccessors
#' @export
setGeneric("symbol", function(object) {
    standardGeneric("symbol")
})

#' @rdname GeneAccessors
#' @export
setMethod("symbol", "Gene", function(object) {
    object@symbol
})

#' Set the Symbol of a Gene
#'
#' Updates the HUGO symbol of a gene.
#'
#' @inheritParams GeneAccessors
#' @return The updated gene object.
#' @examples
#' library(GenomicRanges)
#' gene <- createProteinCodingGene(
#'   id = "ENSG000001",
#'   symbol = "TP53",
#'   name = "Tumor Protein",
#'   description = "Tumor suppressor gene",
#'   structure = GRanges(
#'     seqnames = "chr17",
#'     ranges = IRanges(start = 7668402, end = 7687550),
#'     strand = "+"
#'   ),
#'   protein_id = "P04637",
#'   protein_sequence = "MVLSPADKTNVKAAWG"
#' )
#' symbol(gene) <- "TP53_UPDATED"
#' @rdname GeneAccessors
#' @export
setGeneric("symbol<-", function(object, value) {
    standardGeneric("symbol<-")
})

#' @rdname GeneAccessors
#' @export
setReplaceMethod("symbol", "Gene", function(object, value) {
    object@symbol <- value
    validObject(object)
    object
})

#' Get the Name of a Gene
#'
#' Retrieves the name of a gene.
#'
#' @inheritParams GeneAccessors
#' @return A character string representing the gene name.
#' @examples
#' library(GenomicRanges)
#' gene <- createProteinCodingGene(
#'   id = "ENSG000001",
#'   symbol = "TP53",
#'   name = "Tumor Protein",
#'   description = "Tumor suppressor gene",
#'   structure = GRanges(
#'     seqnames = "chr17",
#'     ranges = IRanges(start = 7668402, end = 7687550),
#'     strand = "+"
#'   ),
#'   protein_id = "P04637",
#'   protein_sequence = "MVLSPADKTNVKAAWG"
#' )
#' name(gene)
#' @rdname GeneAccessors
#' @export
setGeneric("name", function(object) {
    standardGeneric("name")
})

#' @rdname GeneAccessors
#' @export
setMethod("name", "Gene", function(object) {
    object@name
})

#' Set the Name of a Gene
#'
#' Updates the name of a gene.
#'
#' @inheritParams GeneAccessors
#' @return The updated gene object.
#' @examples
#' library(GenomicRanges)
#' gene <- createProteinCodingGene(
#'   id = "ENSG000001",
#'   symbol = "TP53",
#'   name = "Tumor Protein",
#'   description = "Tumor suppressor gene",
#'   structure = GRanges(
#'     seqnames = "chr17",
#'     ranges = IRanges(start = 7668402, end = 7687550),
#'     strand = "+"
#'   ),
#'   protein_id = "P04637",
#'   protein_sequence = "MVLSPADKTNVKAAWG"
#' )
#' name(gene) <- "TP53 Protein"
#' @rdname GeneAccessors
#' @export
setGeneric("name<-", function(object, value) {
    standardGeneric("name<-")
})

#' @rdname GeneAccessors
#' @export
setReplaceMethod("name", "Gene", function(object, value) {
    object@name <- value
    validObject(object)
    object
})
