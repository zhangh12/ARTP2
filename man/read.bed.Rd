\name{read.bed}
\alias{read.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reading data from binary PLINK files}
\description{
Loads genotype data from PLINK format files .bed, .bim, and .fam.
}
\usage{
read.bed(bed, bim, fam, sel.snps = NULL, sel.subs = NULL, encode012 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{
the name of the bed file. 
}
  \item{bim}{
the name of the bim file. For a SNP without a rs number, use any character (including any white space or '.') in the second column of the bim file. 
}
  \item{fam}{
the name of the fam file. 
}
  \item{sel.snps}{
a character vector of SNPs to be extracted from the plink files. The default is \code{NULL}, i.e., all SNPs are extracted. SNPs could be named by its rs number (e.g. rs1234), or by \code{Chr:Pos} (e.g. 13:234567, or C13P234567) if a rs number is not available. All other naming methods for a SNP are not accepted in current version. 
}
  \item{sel.subs}{
an optional character vector specifying a subset of subject IDs to be extracted from the plink files. These IDs should be matched with the second column of \code{fam} files. The default is \code{NULL}, i.e., all subjects are extracted. 
}
  \item{encode012}{
logical. Encoding the genotypes using 0/1/2 if \code{TRUE}, or using symbols of the reference and effect alleles if \code{FALSE}. The default is \code{TRUE}. 
}
}
\value{
A data frame of genotypes of specified subjects in the plink files. For a SNP in \code{sel.snps} specified in the format \code{Chr:Pos}, e.g. 13:234567, it will be named to be C13P234567 in the returned data frame. 
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
