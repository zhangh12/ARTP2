\name{exclude.snps}
\alias{exclude.snps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Excluding SNPs from specified regions}
\description{
\code{exclude.snps} is used to exclude SNPs from specified regions, for example, the GWAS hits and their neighbors. 
}
\usage{
exclude.snps(stat, excluded.regions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{
a data frame containing at least two columns \code{Chr} and \code{Pos}, i.e., the chromosome numbers and base-pair position (bp units) of SNPs. 
}
  \item{excluded.regions}{
a data frame specifying the regions to be excluded. It must contains columns \code{Chr}, \code{Start}, and \code{End}. The unit is base-pair (bp). SNPs within \code{[Start, End]} will be excluded. 
}
}
\value{
A data frame containing a subset of \code{stat} after excluding SNPs in specified regions. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

