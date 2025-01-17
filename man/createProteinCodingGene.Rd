% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_ConstructorF.R
\name{createProteinCodingGene}
\alias{createProteinCodingGene}
\title{Create a Protein Coding Gene Object}
\usage{
createProteinCodingGene(
  id,
  symbol,
  name,
  description,
  structure,
  protein_id,
  protein_sequence
)
}
\arguments{
\item{id}{Character. Gene identifier.}

\item{symbol}{Character. Gene symbol.}

\item{name}{Character. Gene name.}

\item{description}{Character. Gene description.}

\item{structure}{GRanges. Gene structure data.}

\item{protein_id}{Character. Protein identifier.}

\item{protein_sequence}{Character. Protein sequence.}
}
\value{
An instance of the ProteinCodingGene class.
}
\description{
Constructor function for creating a ProteinCodingGene object.
}
\examples{
library(GenomicRanges)
structure <- GRanges(
  seqnames = "chr17",
  ranges = IRanges(start = 7668402, end = 7687550),
  strand = "+"
)
tp53 <- createProteinCodingGene(
  id = "ENSG000001",
  symbol = "TP53",
  name = "Tumor Protein",
  description = "Tumor suppressor gene",
  structure = structure,
  protein_id = "P04637",
  protein_sequence = "MVLSPADKTNVKAAWG"
)
print(tp53)
}
