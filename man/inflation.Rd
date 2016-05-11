\name{inflation}
\alias{inflation}
\title{
Computing the genomic control inflation factor from a given numeric vector
}
\description{
This function gives the inflation of a set of p-values (or 1-df chi-squared statistics). 
}
\usage{
inflation(p, is.p, na.rm = FALSE)
}

\arguments{
  \item{p}{
a numeric vector. It can contain p-values or 1-df chi-squared statistics. 
}
  \item{is.p}{
a logical value indicating whether \code{p} contain p-values or a vector of 1-df chi-squared statistics. 
}
  \item{na.rm}{
a logical value indicating whether \code{NA} values should be stripped before the computation proceeds. }
}
\value{
The inflation factor is returned. 
}
\seealso{
\code{\link{meta}}
}
\examples{
set.seed(0)
p <- runif(10)
inflation(p, TRUE)
}

