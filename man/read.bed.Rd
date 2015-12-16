\name{read.bed}
\alias{read.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data input}
\description{
Loads genotype data from plink format files .bed, .bim, and .fam.
}
\usage{
read.bed(bed, bim, fam, sel.snps = NULL, sel.subs = NULL)
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
a character vector of SNPs to be extracted from the plink files. The default is \code{NULL}, i.e., all SNPs are extracted. 
}
  \item{sel.subs}{
an optional vector specifying a subset of subject IDs to be extracted from the plink files. These IDs should be matched with the second column of \code{fam} files. The default is \code{NULL}, i.e., all subjects are extracted. 
}
}
\value{
A data frame of genotypes of specified subjects in the plink files. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# Load the sample data

bed <- system.file("extdata", package = 'ARTP2', 'chr1.bed')
bim <- system.file("extdata", package = 'ARTP2', 'chr1.bim')
fam <- system.file("extdata", package = 'ARTP2', 'chr1.fam')

## first five SNPs
b <- read.table(bim, header = FALSE, as.is = TRUE, nrows = 5)
## first 50 subjects
f <- read.table(fam, header = FALSE, as.is = TRUE, nrows = 50)
geno <- read.bed(bed, bim, fam, sel.snps = b[, 2], sel.subs = f[, 2])

dim(geno) # 50 x 5


}
