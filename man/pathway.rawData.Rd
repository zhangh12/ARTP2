\name{pathway.rawData}
\alias{pathway.rawData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ARTP test for raw data}
\description{
Calculate gene and pathway p-values using the ARTP test and raw genotype data
}
\usage{
pathway.rawData(formula, data, pathway, family, lambda = 1.0, subset = NULL, options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class \code{\link[stats]{formula}}: a symbolic description of basic risk model to be fitted. Only the outcome and covariates are included. See more details of formula in \code{\link[stats]{glm}}. 
}
  \item{data}{
a data frame containing the variables specified in \code{formula} and genotypes of SNPs. 
}
  \item{pathway}{
a character of the name of file containing definition of a pathway. It must be able to be read by \code{\link[utils]{read.table}} and have columns called "SNP", "Gene", "Chr". It also can be a data frame with the three columns. 
}
  \item{family}{
a character taking values of "\code{gaussian}" or "\code{binomial}". 
}
  \item{lambda}{
a numeric specifying inflation factor. The default is 1.0. 
}
  \item{subset}{
an optional vector specifying a subset of observations in \code{data}. The default is \code{NULL}, i.e., all observations are used. 
}
  \item{options}{
a list of options to control the test procedure. If \code{NULL}, default options will be used. See \code{\link{options}}.
}
}
\details{
This function computes gene and pathway p-values when raw genotype data is available. The ARTP test modified from Yu et al. (2009) and AdaJoint test from Zhang et al. (2014) are released with this package. ARTP is the Adaptive Rank Truncated Product test.
}
\value{
\code{pathway.rawData} returns an object of class \code{"ARTP3"}. It is a list containing the following components: 
\item{pathway.pvalue}{final pathway p-value accounting for multiple comparisons. }
\item{gene.pvalue}{a data frame containing gene name, number of SNPs in the gene that were included in the analysis, chromosome name, and the p-value for the gene accounting for multiple comparisons. }
\item{pathway}{a data frame defining pathway that was actually tested after various filters applied. }
\item{model}{a list containing detailed information of selected SNPs in each gene. }
\item{most.sig.genes}{a character vector of genes selected by \code{ARTP3}. They are the most promising candidates, although its statistical significance is not guaranteed. }
\item{deleted.snps}{a data frame containing SNPs excluded from the analysis and their reasons. }
\item{deleted.genes}{a data frame containing genes excluded from the analysis and their reasons. Set \code{options$rm.gene.subset} to be \code{FALSE} to include all genes. }
\item{options}{a list of options used in the analysis. See \code{\link{options}}. }
\item{accurate}{\code{TRUE} if \code{options$nperm} is large enougth to accurately estimate p-values, i.e., if the criteria 
\code{sqrt(pvalue*(1-pvalue)/nperm)/pvalue < 0.1}
is satisfied.
}
\item{setup}{a list containing necessary input for \code{\link{pathway.warm.start}}. It can be written to a file by using the function \code{\link[base]{save}}, then its path can be the input of \code{\link{pathway.warm.start}}. It also contains a data frame of outcome and covariates that are specified in \code{formula} (\code{setup$yx}), a data frame of genotypes of SNPs in \code{pathway} (\code{setup$raw.geno}), and a formula object \code{setup$formula} corresponding to \code{setup$yx}, if \code{options$keep.geno} is \code{TRUE}. }
}
\references{
Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. 
Genet Epidemiol 33(8): 700 - 709

Zhang H, Shi J, Liang F, Wheeler W, Stolzenberg-Solomon R, Yu K.  (2014) A fast multilocus test with adaptive SNP selection for large-scale genetic association studies. European Journal of Human Genetics, 22, 696 - 702
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{options}}, \code{\link{pathway.warm.start}}, \code{\link{pathway.summaryData}}, \code{\link{raw}}. 
}
\examples{

library(ARTP3)

# Load the sample data
data(raw, package = "ARTP3")
head(raw[, 1:6])

# Define the formula of base risk model
formula <- formula(case_control ~ sex + age + bmi + factor(study))

# Path of a build-in file containing pathway definition
pathway <- system.file('extdata', package = "ARTP3", "pathway.txt.gz")

# binary outcome
family <- "binomial"

# Set the options. 
# Accumulate signal from the top 5 SNPs in each gene
# 1e5 replicates of resampling to estimate the p-value
options <- list(inspect.snp.n = 5, nperm = 1e5, 
                maf = .01, HWE.p = 1e-6, 
                gene.R2 = .9, 
                id.str = "unique-pathway-id", 
                out.dir = getwd(), save.setup = FALSE)

## pathway test
## can take a few minutes
# ret <- pathway.rawData(formula, data = raw, pathway, family, options = options)
# 
# ret$pathway.pvalue 
## [1] 0.02543975 # Mac OS
## [1] 0.02312977 # Windows
# head(ret$gene.pvalue)
# table(ret$deleted.snps$reason)
# head(ret$deleted.genes)


}
