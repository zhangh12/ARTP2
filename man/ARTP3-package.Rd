\name{ARTP3-package}
\alias{ARTP3-package}
\alias{ARTP3}
\docType{package}
\title{
Calculate gene and pathway p-values using the ARTP Test and summary data
}
\description{
Calculate gene and pathway p-values using the ARTP Test and summary data
}
\details{
\tabular{ll}{
Package: \tab ARTP3\cr
Type: \tab Package\cr
Version: \tab 0.8.10\cr
Date: \tab 2015-10-27\cr
License: \tab MIT + file LICENSE\cr
}

It is increasingly recognized that pathway analyses-a joint test of association between the outcome and a group of single nucleotide polymorphisms (SNPs) within a biological pathway could potentially complement single-SNP analysis and provide additional insights for the genetic architecture of complex diseases. Building upon existing P-value combining methods, we propose a class of highly flexible pathway analysis approaches based on an adaptive rank truncated product statistic that can effectively combine evidence of associations over different SNPs and genes within a pathway. The statistical significance of the pathway-level test statistics is evaluated using a highly efficient permutation algorithm that remains computationally feasible irrespective of the size of the pathway and complexity of the underlying test statistics for summarizing SNP- and gene-level associations. 

The main functions in this package are
\code{\link{pathway.summaryData}} when only summary level data are available,
\code{\link{pathway.rawData}} when genotype data are available, and
\code{\link{pathway.warm.start}} for computing gene and pathway p-values when previously save information is available.


}
\author{
Han Zhang, Bill Wheeler, Kai Yu, Yifan Yang

Maintainer: Han Zhang <han.zhang2@nih.gov>
}
\references{
Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. 
Genet Epidemiol 33(8): 700 - 709
}

\keyword{pathway}

