\name{read.bed}
\alias{read.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data input}
\description{
Loads genotype data from plink format files .bed, .bim, and .fam.
}
\usage{
read.bed(bed, bim, fam, sel.snps, sel.subs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{
the name of the bed file. 
}
  \item{bim}{
the name of the bim file. 
}
  \item{fam}{
the name of the fam file. 
}
  \item{sel.snps}{
a character vector of SNPs to be extracted from the plink files. 
}
  \item{sel.subs}{
an optional vector specifying a subset of observations to be extracted from the plink files. The default is \code{NULL}, i.e., all observations are extracted. 
}
}
\value{
A data frame of genotypes of specified observations in the plink files. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# Load the sample data

}
