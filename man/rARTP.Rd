\name{rARTP}
\alias{rARTP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ARTP test for raw data}
\description{
Calculate gene and pathway p-values using the ARTP test and raw genotype data
}
\usage{
rARTP(formula, data, pathway, family, geno.files = NULL, lambda = 1.0, 
      subset = NULL, options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class \code{\link[stats]{formula}}: a symbolic description of basic risk model to be fitted. Only the outcome and covariates are included. See more details of formula in \code{\link[stats]{glm}}. 
}
  \item{data}{
a data frame containing the variables specified in \code{formula}. If \code{geno.files} is \code{NULL}, then it also contains genotypes. 
}
  \item{pathway}{
a character of the name of file containing definition of a pathway. 
It must be able to be read by \code{\link[utils]{read.table}} and have columns called
 \code{SNP}, \code{Gene}, \code{Chr}. 
It also can be a data frame with the three columns. 
The \code{SNP} column can also have values of the form \code{loc1-loc2}, where \code{loc1} and \code{loc2} are
base-pair locations denoting a region of SNPs to use.
}
  \item{family}{
a character taking values of \code{'gaussian'} or \code{'binomial'}. 
}
  \item{geno.files}{
a character vector containing paths of plain text files containing the genotype data. 
Those files can be compressed as gz files and are able to be read by \code{\link[utils]{read.table}}. 
It can be a data frame with columns \code{bed}, \code{bim}, and \code{fam}. 
The data frame contains paths of (multiple sets of) PLINK files containing the genotype data. 
It can be \code{NULL} if all genotype data are already been put in \code{data}. 
}
  \item{lambda}{
a numeric specifying inflation factor. The default is 1.0. 
}
  \item{subset}{
an optional integer vector specifying a subset of observations in \code{data}. The default is \code{NULL}, i.e., all observations are used. 
}
  \item{options}{
a list of options to control the test procedure. If \code{NULL}, default options will be used. See \code{\link{options}}.
}
}
\details{
This function computes gene and pathway p-values when raw genotype data is available. The ARTP test modified from Yu et al. (2009) and AdaJoint test from Zhang et al. (2014) are released with this package. ARTP is the Adaptive Rank Truncated Product test.

The raw (i.e. individual-level) genotype data, can be encoded as 0, 1, or 2 (counts of effect allele), or any quantitative values (e.g., output from genotype imputation program). 
}
\value{
\code{rARTP} returns an object of class \code{ARTP2}. It is a list containing the following components: 
\item{pathway.pvalue}{final pathway p-value accounting for multiple comparisons. }
\item{gene.pvalue}{a data frame containing gene name, number of SNPs in the gene that were included in the analysis, chromosome name, and the p-value for the gene accounting for multiple comparisons. }
\item{pathway}{a data frame defining the pathway that was actually tested after various filters applied. }
\item{model}{a list containing detailed information of selected SNPs in each gene. }
\item{most.sig.genes}{a character vector of genes selected by \code{ARTP2}. They are the most promising candidates, although their statistical significance is not guaranteed. }
\item{deleted.snps}{a data frame containing SNPs excluded from the analysis and their reasons. }
\item{deleted.genes}{a data frame containing genes excluded from the analysis because they are subsets of other remaining genes. Set \code{options$rm.gene.subset} to be \code{FALSE} to include all genes even if they are subsets of other genes. }
\item{options}{a list of options used in the analysis. See \code{\link{options}}. }
\item{accurate}{\code{TRUE} if \code{options$nperm} is large enougth to accurately estimate p-values, i.e., if the criteria 
\code{sqrt(pvalue*(1-pvalue)/nperm)/pvalue < 0.1}
is satisfied.
}
\item{setup}{a list containing necessary input for \code{\link{warm.start}}. It can be written to a file by using the function \code{\link[base]{save}}, then its path can be the input of \code{\link{warm.start}}. It also contains a data frame of outcome and covariates that are specified in \code{formula} (\code{setup$yx}), a data frame of genotypes of SNPs in \code{pathway} (\code{setup$raw.geno}), and a formula object \code{setup$formula} corresponding to \code{setup$yx}, if \code{options$keep.geno} is \code{TRUE}. }
}
\references{
Zhang H, Wheeler W, Hyland LP, Yang Y, Shi J, Chatterjee N, Yu K. (2016) A powerful procedure for pathway-based meta-analysis using summary statistics identifies 43 pathways associated with type II diabetes in European populations. PLoS Genetics 12(6): e1006122

Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. Genet Epidemiol 33(8): 700 - 709

Zhang H, Shi J, Liang F, Wheeler W, Stolzenberg-Solomon R, Yu K.  (2014) A fast multilocus test with adaptive SNP selection for large-scale genetic association studies. European Journal of Human Genetics 22: 696 - 702
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{options}}, \code{\link{warm.start}}, \code{\link{sARTP}}, \code{\link{example}}. 
}
\examples{

library(ARTP2)

## Load the sample data
data(data, package = "ARTP2")
head(data[, 1:7])

## Load a build-in data frame containing pathway definition
## it can also be the path of the file
data(pathway, package = "ARTP2")
head(pathway)

## Define the formula of base risk model
formula <- formula(case_control ~ sex + age + bmi + factor(study))

## binary outcome
family <- "binomial"

## Set the options. 
## Accumulate signal from the top 5 SNPs in each gene
## 1e5 replicates of resampling to estimate the p-value
options <- list(inspect.snp.n = 5, nperm = 1e5, 
                maf = .01, HWE.p = 1e-6, 
                gene.R2 = .9, 
                id.str = "unique-pathway-id", 
                out.dir = getwd(), save.setup = FALSE)

## pathway test, can take a while
## data contains outcome, covariates and genotypes
# ret1 <- rARTP(formula, data = data, pathway, family, options = options)

# ret1$pathway.pvalue 
## [1] 0.03218968 # Mac OS
## [1] 0.02188978 # Linux with 1 thread
## [1] 0.03455965 # Linux with 32 threads

## Mac OS
# head(ret1$gene.pvalue)
##     Gene Chr N.SNP      Pvalue
## 1  USP30  12    18 0.001319987
## 2  DCAF7  17     9 0.071644284
## 3   CANX   5    13 0.266337337
## 4  SOX12  20    15 0.349406506
## 5 CDKN2C   1     6 0.358031420
## 6   FEN1  11     4 0.415345847

## Linux with 1 thread
# head(ret1$gene.pvalue)
##     Gene Chr N.SNP      Pvalue
## 1  USP30  12    18 0.000899991
## 2  DCAF7  17     9 0.070219298
## 3   CANX   5    13 0.269772302
## 4  SOX12  20    15 0.350061499
## 5 CDKN2C   1     6 0.357766422
## 6   FEN1  11     4 0.414760852

## Linux with 32 threads
# head(ret1$gene.pvalue)
##     Gene Chr N.SNP      Pvalue
## 1  USP30  12    18 0.001454985
## 2  DCAF7  17     9 0.070379296
## 3   CANX   5    13 0.266927331
## 4  SOX12  20    15 0.350481495
## 5 CDKN2C   1     6 0.357701423
## 6   FEN1  11     4 0.414425856

# table(ret1$deleted.snps$reason)
# head(ret1$deleted.genes)


##################################################
## Another way to use this function
## Load a vector 'geno' containing file names of genotype
data(geno, package = 'ARTP2')

## Set the paths of genotype files
## in this example, each file contains SNPs in a gene
geno.files <- system.file("extdata", package = "ARTP2", geno)

## data contains outcome, covariates
## Genotypes are instead included in files specified in geno.files
## geno.files are plain text files (or .gz file), which can be read by read.table
# ret2 <- rARTP(formula, data = data[, 2:6], pathway, family, geno.files, 
#               options = options)
# ret2$pathway.pvalue == ret1$pathway.pvalue


##################################################
## The third way
## Genotypes are instead stored as binary PLINK files (bed, bim, and fam)
bed <- system.file("extdata", package = "ARTP2", "raw.bed")
bim <- system.file("extdata", package = "ARTP2", "raw.bim")
fam <- system.file("extdata", package = "ARTP2", "raw.fam")
geno.files <- data.frame(fam, bim, bed, stringsAsFactors = FALSE)

## a column SUBID must be included in data, in this example, first column is SUBID
# ret3 <- rARTP(formula, data = data[, 1:6], pathway, family, geno.files, 
#               options = options)
# ret3$pathway.pvalue == ret1$pathway.pvalue


}
