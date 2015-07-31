\name{options.default}
\alias{options.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
options.default
}
\description{
A function to return default \code{\link{options}} in \code{\link{pathway.summaryData}} and \code{\link{pathway.rawData}}. This page is for illustration. Users do not need to call this fuction explicitly.
}
\usage{
options.default()
}
\value{
A list of \code{\link{options}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{options}}
}
\examples{
options <- options.default()
str(options)
names(options)
}
