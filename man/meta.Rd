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
A data frame of summary statistics from meta-analysis, and of individual studies specified in \code{summary.files}. 
}
\references{
  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    
  }
\examples{
  
  
}

