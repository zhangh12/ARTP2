\name{pathway.warm.start}
\alias{pathway.warm.start}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pathway.warm.start
}
\description{
This function is designed to accelerate the \code{ARTP3} test in practice. It uses pre-calculated and reusable statistics as input and allow the users to try different testing configuration more efficiently. See \code{Details} for more information. 
}
\usage{
pathway.warm.start(setup.file, nperm = NULL, inflation.factor = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setup.file}{character string of file name saving the setup created by \code{\link{pathway.summaryData}} or \code{\link{pathway.rawData}}. }
  \item{nperm}{the number of permutations. If it is \code{NULL} as default, then the value in the saved setup will be used. See \code{Details}. }
  \item{inflation.factor}{inflation factor lambda to be adjusted in pathway analysis. The default is 1.0 (NULL). }
}
\details{
An \code{ARTP3} test has two major steps in testing an association. The first step applys data clean criteria and creates necessary and reusable statistics, which can be time-consuming for large pathways. The second step performs the testing procedure to estimate the pathway or gene-level p-value. \code{pathway.warm.start} focuses on the second step. 

The first step can be done by using \code{\link{pathway.summaryData}} or \code{\link{pathway.rawData}} if their \code{options$only.setup} is set as \code{TRUE}. A list \code{setup} in their output can be \code{\link[base]{save}} to a local file, and its path can be used as the input of \code{pathway.warm.start}. With \code{pathway.warm.start}, users can try different configurations to perform a various of tests allowed by the \code{ARTP3} framework, but avoid long waiting time for data cleaning. Commonly used options in \code{setup$options} include \code{method}, \code{inspect.snp.n}, \code{inspect.gene.n}, \code{nperm}, \code{inflation.factor}, etc. 

Note that both \code{\link{pathway.summaryData}} and \code{\link{pathway.rawData}} can produce the final p-value directly if \code{options$only.setup} is \code{FALSE}. 

The \code{setup.file} is supposed to be the path of a \code{".rda"} file containing a list called \code{setup}. If \code{nperm} is \code{NULL}, then it will be set as \code{setup$options$nperm}. If \code{inflation.factor} is \code{NULL}, then it will be set as 1.0 as default. These two options are the most useful ones in using \code{pathway.warm.start} so we highlight them in the interface. Users can modify any option values in \code{setup$options} directly to get more controls of the testing procedure. See \code{\link{options}} for more details about setting \code{setup$options}. 

Except for \code{setup$options}, all other components in \code{setup} cannot be modified by users. Modified \code{setup} should be saved to local file \code{setup.file} by \code{\link[base]{save}} so that \code{pathway.warm.start} can make use of the most updated configuration. 
}
\value{
\code{pathway.warm.start} returns an object of class \code{"ARTP3"}. It is a list containing the following components: 
\item{pathway.pvalue}{final pathway p-value accounting for multiple comparisons.}
\item{gene.pvalue}{a data frame containing gene name, number of SNPs in the gene that were included in the analysis, chromosome name, and the p-value for the gene accounting for multiple comparisons.}
\item{pathway}{a data frame defining pathway that was actually tested after various filters applied. }
\item{model}{a list containing detailed information of selected SNPs in each gene. }
\item{most.sig.genes}{a character vector of genes selected by \code{ARTP3}. They are the most promising candidates, although its statistical significance is not guaranteed. }
\item{deleted.snps}{a data frame containing SNPs excluded from the analysis and their reasons.}
\item{options}{a list of options used in the analysis. See \code{\link{options}}}
\item{test.timing}{timing information (in sec)}
\item{accurate}{\code{TRUE} if \code{options$nperm} is large enougth to accurately estimate p-values, i.e., if the criteria 
\code{sqrt(pvalue*(1-pvalue)/nperm)/pvalue < 0.1}
is satisfied.
}
}
\references{
Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N. (2009) Pathway analysis by adaptive combination of P-values. 
Genet Epidemiol 33(8): 700 - 709
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pathway.summaryData}}, \code{\link{pathway.rawData}}, \code{\link{options}}
}
\examples{

## firstly, run the example in pathway.rawData and get an ARTP3 object "ret"

# setup <- ret$setup
## modify or specify the method
# setup$options$method <- 2
# setup$options$inspect.snp.n <- 5
# setup.file <- "path/to/your/file.rda"
# save(setup, file = setup.file)

# ret1 <- pathway.warm.start(setup.file, nperm = 1e5)
# ret2 <- pathway.warm.start(setup.file, nperm = 1e4, inflation.factor = 1.08)


}

