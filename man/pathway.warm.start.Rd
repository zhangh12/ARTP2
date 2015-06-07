\name{pathway.warm.start}
\alias{pathway.warm.start}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pathway.warm.start
}
\description{
Calculate gene and pathway p-values using the ARTP test. The input is created from summary data by calling \code{\link{pathway.summaryData}} with \code{options$only.setup = TRUE}. 
}
\usage{
pathway.warm.start(setup.file, nperm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setup.file}{character string of file name saving the setup created by \code{\link{pathway.summaryData}}. }
  \item{nperm}{the number of permutations. If it is \code{NULL} as default, then the value in the saved setup will be used. See \code{Details}. }
}
\details{
The setup.file is supposed to be a \code{".rda"} R file containing a list called \code{setup}. If \code{nperm} is \code{NULL}, then it will be set as \code{setup$options$nperm}. The user can modify the option values in \code{setup$options} to get more control of the testing procedure. 
}
\value{
\code{pathway.warm.start} returns an object of class \code{"ARTP3"}. It is a list containing the following components: 
\item{pathway.pvalue}{final pathway p-value accounting for multiple comparisons.}
\item{gene.pvalue}{a data frame containing gene name, number of SNPs in the gene that were included in the analysis, chromosome name, and the p-value for the gene accounting for multiple comparisons.}
\item{pathway}{a data frame defining pathway that was actually tested after various filters applied. }
\item{deleted.snps}{a data frame containing SNPs excluded from the analysis and their reasons.}
\item{options}{a list of options used in the analysis. See \code{\link{options}}}
\item{timestamp}{a vector of timing information.}
\item{accurate}{\code{TRUE} if \code{options$nperm} is large enougth to accurately estimate p-values, i.e., if the criteria 
\code{sqrt(pvalue*(1-pvalue)/nperm)/pvalue < 0.1}
is satisfied.
}
}
\references{
Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. 
Genet Epidemiol 33(8): 700 - 709
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pathway.summaryData}}
}
\examples{
## example

}

