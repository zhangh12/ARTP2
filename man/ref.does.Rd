\name{ref.does}
\alias{ref.does}
\docType{data}
\title{
  Example reference genotype data coded as expected dosages. 
}
\description{
  A list with names \code{ref.geno} and \code{allele.info} containing genotypes
 coded as expected dosages and the corresponding allele information for each SNP.
}
\usage{data(ref.does)}

\examples{
  data(ref.does)
  head(ref.does$ref.geno[, 1:5])
  head(ref.does$allele.info[, 1:5])
}
\keyword{datasets}
