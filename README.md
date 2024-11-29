# GeneRepresentationPackage

The `GeneRepresentationPackage` provides tools for creating and working with S4 classes 
representing genes, including protein-coding genes, long non-coding RNAs (lncRNAs) 
and microRNAs (miRNAs). These classes inherit from a virtual `Gene` class.

## Features
- S4 class structure for genes and their subtypes
- Constructor functions for creating gene objects
- Accessor methods to retrieve and modify attributes
- Functionality to compute product lengths for different gene types

## Installation
Clone this repository and use the following command to install the package:

```r
devtools::install_github("lorenzomonticelli/GeneRepresentationPackage")
```