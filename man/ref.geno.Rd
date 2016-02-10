\name{ref.geno}
\alias{ref.geno}
\docType{data}
\title{
  A dataset used in example of \code{pathway.summaryData}. 
}
\description{
  A data frame \code{ref.geno} with genotypes of 503 observations on 2654 SNPs. The genotypes are encoded by SNPs' minor and major alleles so that this data frame can be used as \code{reference} of function \code{pathway.summaryData}. 
}
\usage{data(ref.geno)}

\examples{
  data(ref.geno)
  head(ref.geno[, 1:5])
}
\keyword{datasets}
