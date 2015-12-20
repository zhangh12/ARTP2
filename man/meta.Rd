\name{meta}
\alias{meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  meta
}
\description{
  Perform the single-marker meta-analysis. 
}
\usage{
  meta(summary.files, lambda = NULL, sel.snps = NULL, only.meta = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summary.files}{
a character vector of file names containing the summary results of SNPs included in one or multiple studies. Each file must be able to be read by \code{\link[utils]{read.table}}. Each file must have columns called "SNP", "RefAllele", "EffectAllele", "BETA", and at least one of "SE", "P".
}
  \item{lambda}{
a numeric vector of inflation factors. Each file in \code{summary.files} should have one inflation factor specified in \code{lambda}. \code{NULL} if inflation is not adjusted. 
}
  \item{sel.snps}{
a character vector of SNPs to be used in meta-analysis. The default is \code{NULL}, i.e., all SNPs are used. 
}
  \item{only.meta}{
\code{TRUE} if do not returned individual summary data. The default is \code{TRUE}. 
}
}
\details{
  The inverse-variance method is used by assuming a fixed effect model. The standard error is rescaled by \code{sqrt{lambda}}. 
}
\value{
\code{meta} returns a list containing 
\item{meta.stat}{a data frame of summary statistics from meta-analysis, and of individual studies specified in \code{summary.files} if \code{only.meta} is \code{FALSE}. }
\item{conf.snps}{a character vector of SNPs with conflictive allele information. }
}
\references{
  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    
  }
\examples{
  
  study1 <- system.file("extdata", package = "ARTP2", "study1.txt.gz")
  study2 <- system.file("extdata", package = "ARTP2", "study2.txt.gz")
  snps <- c('rs13266821', 'rs4824130', 'rs1000047', 'rs1000017', 'rs6066771', 'rs12508128')
  
  m1 <- meta(summary.files = c(study1, study2), lambda = c(1.10, 1.08), sel.snps = snps, only.meta = FALSE)
  m2 <- meta(summary.files = c(study1, study2), lambda = c(1.10, 1.08), sel.snps = snps, only.meta = TRUE)
  
  m1$meta.stat
  m2$meta.stat
  
  is.null(m1$conf.snps)
  
}

