\name{ARTP2-package}
\alias{ARTP2-package}
\alias{ARTP2}
\docType{package}
\title{
Pathway and Gene-Level Association Test
}
\description{
Pathway and gene level association test using raw data or summary statistics. 
}
\details{
\tabular{ll}{
Package: \tab ARTP2\cr
Type: \tab Package\cr
Version: \tab 0.9.22\cr
Date: \tab 2016-02-11\cr
License: \tab GPL-2 | GPL-3\cr
}

It is increasingly recognized that pathway analyses, a joint test of association between the outcome and a group of single nucleotide polymorphisms (SNPs) within a biological pathway, could potentially complement single-SNP analysis and provide additional insights for the genetic architecture of complex diseases. Building upon existing P-value combining methods, we propose a class of highly flexible pathway analysis approaches based on an adaptive rank truncated product statistic that can effectively combine evidence of associations over different SNPs and genes within a pathway. The statistical significance of the pathway-level test statistics is evaluated using a highly efficient permutation algorithm that remains computationally feasible irrespective of the size of the pathway and complexity of the underlying test statistics for summarizing SNP- and gene-level associations. 

The main functions in this package are
\code{\link{sARTP}} when only summary level data are available,
\code{\link{rARTP}} when genotype data are available, and
\code{\link{warm.start}} for computing gene and pathway p-values when previously save information is available.


}
\author{
Han Zhang, Bill Wheeler, Kai Yu, Yifan Yang

Maintainer: Han Zhang <han.zhang2@nih.gov>
}
\references{
Zhang H, Wheeler W, Hyland LP, Yang Y, Shi J, Chatterjee N, Yu K. (2016) A powerful procedure for pathway-based meta-analysis using summary statistics identifies 43 pathways associated with type II diabetes in European populations. PLoS Genetics 12(6): e1006122

Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. 
Genet Epidemiol 33(8): 700 - 709

Zhang H, Shi J, Liang F, Wheeler W, Stolzenberg-Solomon R, Yu K.  (2014) A fast multilocus test with adaptive SNP selection for large-scale genetic association studies. European Journal of Human Genetics 22: 696 - 702

}

\keyword{pathway}

