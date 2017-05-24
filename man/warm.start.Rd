\name{warm.start}
\alias{warm.start}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing ARTP2 test with warm start
}
\description{
This function is designed to accelerate the \code{ARTP2} test in practice. It uses pre-calculated and reusable statistics as input and allow the users to try different testing configuration more efficiently. See \code{Details} for more information. 
}
\usage{
warm.start(setup, nperm = NULL, lambda = 1.0, nthread = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setup}{an R object created by \code{\link{sARTP}} or \code{\link{rARTP}}. It is a list containing necessary statistics for computing p-values. }
  \item{nperm}{the number of permutations. If it is \code{NULL} as default, then the value in the saved setup will be used. See \code{Details}. }
  \item{lambda}{inflation factor to be adjusted in pathway analysis. \code{lambda} in this function can only be a single numeric number, which is different from the one in \code{sARTP}. The default is 1.0. }
  \item{nthread}{number of threads to be used in permutation. \code{NULL} if \code{setup$options$nthread} is used. }
}
\details{
An \code{ARTP2} test has two major steps in testing an association. The first step applys data clean criteria and creates necessary and reusable statistics, which can be time-consuming for large pathways. The second step performs the testing procedure to estimate the pathway or gene-level p-value. \code{warm.start} focuses on the second step. 

The first step can be done by using \code{\link{sARTP}} or \code{\link{rARTP}} if their \code{options$only.setup} is set as \code{TRUE}. Their output object, \code{setup}, can be used as the first argument of \code{warm.start}. With \code{warm.start}, users can try different configurations to perform various tests allowed by the \code{ARTP2} framework, but avoid long waiting time for data cleaning. Commonly used options in \code{setup$options} include \code{method}, \code{inspect.snp.n}, \code{inspect.gene.n}, \code{nperm}, etc. 

Note that both \code{\link{sARTP}} and \code{\link{rARTP}} can produce the final p-value directly if \code{options$only.setup} is \code{FALSE}. 

The \code{setup} is supposed to have all components defined in \code{\link{sARTP}} and \code{\link{rARTP}}. If \code{nperm} is \code{NULL}, then it will be set as \code{setup$options$nperm}. The users can also pass \code{lambda} if a second round genomic control is needed. However, unlike in \code{\link{sARTP}}, \code{lambda} here can only be a single numberic number rather than a vector. Options \code{nperm} and \code{lambda} are the most useful ones in using \code{warm.start} so we highlight them in the interface. Users can modify any option values in \code{setup$options} directly to get more controls of the testing procedure. See \code{\link{options}} for more details about how to set \code{setup$options}. 

Except for \code{setup$options}, all other components in \code{setup} should not be modified by users. 
}
\value{
\code{warm.start} returns an object of class \code{ARTP2}. It is a list containing the following components: 
\item{pathway.pvalue}{final pathway p-value accounting for multiple comparisons.}
\item{gene.pvalue}{a data frame containing gene name, number of SNPs in the gene that were included in the analysis, chromosome name, and the p-value for the gene accounting for multiple comparisons.}
\item{pathway}{a data frame defining the pathway that was actually tested after various filters applied. }
\item{model}{a list containing detailed information of selected SNPs in each gene. }
\item{most.sig.genes}{a character vector of genes selected by \code{ARTP2}. They are the most promising candidates, although their statistical significance is not guaranteed. }
\item{deleted.snps}{a data frame containing SNPs excluded from the analysis and their reasons.}
\item{deleted.genes}{a data frame containing genes excluded from the analysis because they are subsets of other remaining genes. }
\item{options}{a list of options used in the analysis. See \code{\link{options}}}
\item{test.timing}{timing information (in sec)}
\item{accurate}{\code{TRUE} if \code{options$nperm} is large enougth to accurately estimate p-values, i.e., if the criteria 
\code{sqrt(pvalue*(1-pvalue)/nperm)/pvalue < 0.1}
is satisfied.
}
}
\references{
Zhang H, Wheeler W, Hyland LP, Yang Y, Shi J, Chatterjee N, Yu K. (2016) A powerful procedure for pathway-based meta-analysis using summary statistics identifies 43 pathways associated with type II diabetes in European populations. PLoS Genetics 12(6): e1006122

Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. Genet Epidemiol 33(8): 700 - 709

Zhang H, Shi J, Liang F, Wheeler W, Stolzenberg-Solomon R, Yu K.  (2014) A fast multilocus test with adaptive SNP selection for large-scale genetic association studies. European Journal of Human Genetics: 22, 696 - 702

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sARTP}}, \code{\link{rARTP}}, \code{\link{options}}
}
\examples{

## firstly, run the example in sARTP
## users can adjust the second round inflation in warm.start
## the first round inflation can be study-specific and adjusted in rARTP 
## or sARTP

library(ARTP2)
study1 <- system.file("extdata", package = "ARTP2", "study1.txt.gz")
study2 <- system.file("extdata", package = "ARTP2", "study2.txt.gz")
pathway <- system.file("extdata", package = "ARTP2", "pathway.txt.gz")
chr <- 1:22
nchr <- length(chr)
fam <- vector("character", nchr)
bim <- vector("character", nchr)
bed <- vector("character", nchr)
for(i in 1:nchr){
  fam[i] <- system.file("extdata", package = "ARTP2", paste("chr", chr[i], ".fam", sep = ""))
  bim[i] <- system.file("extdata", package = "ARTP2", paste("chr", chr[i], ".bim", sep = ""))
  bed[i] <- system.file("extdata", package = "ARTP2", paste("chr", chr[i], ".bed", sep = ""))
}
reference <- data.frame(fam, bim, bed)
options <- list(inspect.snp.n = 2, nperm = 1e4, 
                maf = .01, HWE.p = 1e-6, 
                gene.R2 = .9, 
                id.str = "unique-pathway-id", 
                out.dir = getwd(), save.setup = FALSE)
                
## different inflation factors are adjusted in two studies
## first round adjustment
lambda <- c(1.10, 1.08)
ncases <- list()
ncontrols <- list()
ncases[[1]] <- c(9580, 2591)
ncontrols[[1]] <- c(53810, 3052)
ncases[[2]] <- 7638
ncontrols[[2]] <- 54319

family <- 'binomial'

## do not run permutation
options$only.setup <- TRUE
## the first round study-specific inflation is adjusted as lambda = c(1.10, 1.08)
# setup <- sARTP(summary.files = c(study1, study2), pathway, family, reference, lambda, 
#                ncases, ncontrols, options = options)

## the two rounds of inflation is adjusted as lambda2 = c(1.17370, 1.15236)
lambda2 <- lambda * 1.067
## run permutation to calculate p-value
options$only.setup <- FALSE
# ret1 <- sARTP(summary.files = c(study1, study2), pathway, family, reference, lambda2, 
#               ncases, ncontrols, options = options)

## or adjust the second round of inflation in warm.start
# ret2 <- warm.start(setup, lambda = 1.067)

# two ways of inflation adjustment should give same results
# ret1$pathway.pvalue == ret2$pathway.pvalue

###############################################################
###############################################################
## modify or specify the method
# setup$options$method <- 2
# setup$options$inspect.snp.n <- 3

## nthread = 2 for Linux only
## nthread will be reset to 1 under Windows and Mac OS
# ret3 <- warm.start(setup, nperm = 1e5, nthread = 2)

}

