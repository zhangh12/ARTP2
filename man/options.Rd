\name{options}
\alias{options}
\docType{data}
\title{options}
\description{
The list to describe the options that are used in \code{\link{pathway.summaryData}} and \code{\link{pathway.rawData}}. It will be set by function \code{\link{options.default}} by default. 
}
\format{
  The format is a list. 
  
  \bold{Options for testing an association: }
  \describe{
  \item{\code{out.dir}}{output directory for temporary files. The default is the working directory \code{\link[base]{getwd}}. }
  \item{\code{id.str}}{character string that is appended to temporary file names. The default is "PID". }
  \item{\code{method}}{1 = AdaJoint, 2 = AdaJoint2, 3 = ARTP. The default is 3. }
  \item{\code{nperm}}{the number of permutations. The default is 1E5. }
  \item{\code{inflation.factor}}{inflation factor lambda to be adjusted in pathway analysis. The default is 1.0. }
  \item{\code{nthread}}{the number of threads for multi-threaded processors in Unix/Linux OS. The default is \code{detectCores()} to use all available processors. }
  \item{\code{seed}}{integer for random number generation. The default is 1. }
  }
  
  \bold{Options for controlling data cleaning: }
  \describe{
    \item{\code{out.dir}}{output directory for temporary files. The default is the working directory \code{\link[base]{getwd}}. }
    \item{\code{id.str}}{character string that is appended to temporary file names. The default is "PID". }
    \item{\code{snp.miss.rate}}{any SNP with missing rate greater than \code{snp.miss.rate} will be removed from the analysis. The default is 0.05. }
    \item{\code{maf}}{any SNP with minor allele frequency less than \code{maf} will be removed from the analysis. The default is 0.05. }
    \item{\code{HWE.p}}{any SNP with HWE exact p-value less than \code{HWE.p} will be removed from the analysis. The test is applied to the reference data. The default is 1E-5. }
    \item{\code{gene.R2}}{a number between 0 and 1 to filter out SNPs that are highly correlated within each gene. The \code{\link[stats]{cor}} function will be called to compute the R^2 values between each pair of SNPs and remove one SNP with lower MAF in each pair with R^2 greater than \code{gene.R2}. The default is (1) \code{chr.R2} if only \code{chr.R2} is set by the users; (2) 0.95 if \code{chr.R2} is not set by the users. }
    \item{\code{chr.R2}}{a number between 0 and 1 to filter out SNPs that are highly correlated within each chromosome. The \code{\link[stats]{cor}} function will be called to compute the R^2 values between each pair of SNPs and remove one SNP with lower MAF in each pair with R^2 greater than \code{chr.R2}. The default is (1) \code{gene.R2} if only \code{gene.R2} is set by the users; (2) 0.95 if \code{gene.R2} is not set by the users. }
    \item{\code{gene.miss.rate}}{threshold to remove genes based on their missing rate. Genes with missing rate greater than \code{gene.miss.rate} will be removed from the analysis. The missing rate is calculated as the number of subjects with at least one missing genotype among all SNPs in the gene divided by the total number of subjects. The default is 1.0. }
    \item{\code{rm.gene.subset}}{\code{TRUE} to remove genes which are subsets of other genes. The default is \code{TRUE}. }
    \item{\code{delete}}{\code{TRUE} to delete temporary files containing the test statistics for each gene. The default is \code{TRUE}. }
    \item{\code{print}}{\code{TRUE} to print information to the console. The default is \code{TRUE}. }
    \item{\code{tidy}}{the data frame \code{deleted.snps} in the returned object of \code{\link{pathway.summaryData}} containing information of SNPs excluded from the analysis and their reasons. Possible reason codes include \code{"RM_BY_USER"}, \code{"NO_SUM_STAT"}, \code{"NO_REF"}, \code{"SNP_MISS_RATE"}, \code{"SNP_LOW_MAF"}, \code{"SNP_CONST"}, \code{"SNP_HWE"}, \code{"GENE_R2"}, \code{"HUGE_GENE_R2"}, \code{"CHR_R2"}, \code{"GENE_MISS_RATE"}, \code{"GENE_SUBSET"}, \code{"CONF_ALLELE_INFO"}. Set \code{tidy} as \code{TRUE} to hide the SNPs with codes \code{"NO_SUM_STAT"} and \code{"NO_REF"}. The default is \code{TRUE}. }
    \item{\code{save.setup}}{\code{TRUE} to save necessary data, e.g., working options, observed scores and covariance matrix, to local to repeat the analysis more quicly (skip loading and filtering data). It will be set to be \code{TRUE} if \code{only.setup} is \code{TRUE}. The default is \code{TRUE}. }
    \item{\code{path.setup}}{character string of file name to save the setup for \code{\link{pathway.warm.start}} if \code{save.setup} is \code{TRUE}. The default is \code{NULL} so that it is set as \code{paste(out.dir, "/setup.", id.str, ".rda", sep = "")}. }
    \item{\code{only.setup}}{\code{TRUE} if only the setup is needed while the testing procedure is not. The R code to create the setup uses single thread but the testing procedure can be multi-threaded. The best practice to use \code{ARTP3} on a multi-threaded cluster is to firstly create the setup in single-thread mode, and then call the \code{\link{pathway.warm.start}} to compute the p-values in multiple-thread mode, which uses the saved setup at \code{path.setup} as input. \code{save.setup} will be set to be \code{TRUE} if \code{only.setup} is \code{TRUE}. The default is \code{FALSE}. }
    \item{\code{keep.geno}}{\code{TRUE} if the reference genotypes of SNPs in pathway is returned. The default is \code{FALSE}. }
    \item{\code{seed}}{integer for random number generation. The default is 1. }
    \item{\code{excluded.snps}}{character vector of SNPs to be excluded in the analysis. \code{NULL} if no SNP is always excluded. The default is \code{NULL}. }
    \item{\code{selected.snps}}{character vector of SNPs to be selected in the analysis. \code{NULL} if no SNP is always selected but other fiters may applied. The default is \code{NULL}. }
    \item{\code{excluded.subs}}{character vector of subject IDs to be excluded in the analysis. These IDs must match with those in \code{fam} files in \code{reference}. The default is \code{NULL}. }
    \item{\code{selected.subs}}{character vector of subject IDs to be selected in the analysis. These IDs must match with those in \code{fam} files in \code{reference}. The default is \code{NULL}. }
    \item{\code{excluded.genes}}{character vector of genes to be excluded in the analysis. \code{NULL} if no gene is always excluded. The default is \code{NULL}. }
    }
  
  \bold{Options for handling huge pathways: }
  \describe{
    \item{\code{trim.huge.chr}}{oversized chromosomes could be further trimmed to accelerate the testing procedure. If \code{TRUE} the additional options below are in effect. The default is \code{TRUE}. }
    \item{\code{huge.chr.size}}{a chromosome with number of SNPs larger than \code{huge.chr.size} will be further trimmed if \code{trim.huge.pathway} is \code{TRUE}. The default is 2000. }
    \item{\code{huge.gene.R2}}{more stringent R^2 threshold to filter out SNPs in a gene. Similar to \code{gene.R2}. The default is \code{gene.R2} - 0.05. }
    \item{\code{huge.chr.R2}}{more stringent R^2 threshold to filter out SNPs in a chromosome. Similar to \code{chr.R2}. The default is \code{chr.R2} - 0.05. }
  }
  
  
  \bold{Options for gene-based test: } 
  \describe{
    \item{\code{inspect.snp.n}}{the number of candidate truncation points to inspect the top SNPs in a gene. The default is 5. }
    \item{\code{inspect.snp.percent}}{a value \code{x} between 0 and 1 such that a truncation point will be defined at every \code{x} percent of the top SNPs. The default is 0 so that the truncation points will be 1:\code{inspect.snp.n}. }
  }
  
  \bold{Options for pathway-based test: } 
  \describe{
    \item{\code{inspect.gene.n}}{the number of candidate truncation points to inspect the top genes in the pathway. The default is 10. }
    \item{\code{inspect.gene.percent}}{a value \code{x} between 0 and 1 such that a truncation point will be defined at every \code{x} percent of the top genes. If 0 then the truncation points will be 1:\code{inspect.gene.n}. The default is 0.05. }
  }
}

\details{ 
  \bold{Order of removing SNPs, genes and subjects:} \cr
   1. Apply the option \code{excluded.snps} and \code{selected.snps} if non-NULL. Code: \code{"RM_BY_USER"}.\cr
   2. Remove SNPs without summary statistics in \code{summary.files}. Code: \code{"NO_SUM_STAT"}. \cr
   3. Remove SNPs not in \code{bim} files in \code{reference}. Code: \code{"NO_REF"}. \cr
   4. Remove SNPs with conflictive allele information in summary and reference data. Code: \code{"CONF_ALLELE_INFO"}. \cr
   5. Remove SNPs with high missing rate. Code: \code{"SNP_MISS_RATE"}. \cr
   6. Remove SNPs with low MAF. Code: \code{"SNP_LOW_MAF"}. \cr
   7. Remove constant SNPs. Code: \code{"SNP_CONST"}. \cr
   8. Remove SNPs fail to pass HWE test. Code: \code{"SNP_HWE"}. \cr
   9. Remove highly correlated SNPs within each gene. Code: \code{"GENE_R2"} or \code{"HUGE_GENE_R2"}. \cr
   10. Remove highly correlated SNPs within each chromosome. Code: \code{"CHR_R2"}. \cr
   11. Remove genes with high missing rate. Code: \code{"GENE_MISS_RATE"}. \cr
   12. Remove genes which are subsets of other genes. Code: \code{"GENE_SUBSET"}. \cr
   
  \bold{Example truncation points defined by \code{inspect.snp.n} and \code{inspect.snp.percent}: } 
  Assume the number of SNPs in a gene is 100. Below are examples of the truncation points for different
 values of \code{inspect.snp.n} and \code{inspect.snp.percent}. 

\tabular{rrr}{
inspect.snp.n \tab inspect.snp.percent \tab truncation points \cr
1 \tab 0 \tab 1 \cr
1 \tab 0.05 \tab 5 \cr
1 \tab 0.25 \tab 25 \cr
1 \tab 1 \tab 100 \cr
2 \tab 0 \tab 1, 2 \cr
2 \tab 0.05 \tab 5, 10 \cr
2 \tab 0.25 \tab 25, 50 \cr
2 \tab 1 \tab 100 \cr
3 \tab 0.2 \tab 20, 40, 60 
} % END: tabular

}

\seealso{
\code{\link{options.default}}
}

\examples{
options <- options.default()
str(options)
names(options)
}

