
meta <- function(summary.files, lambda = NULL, sel.snps = NULL){
  
  validate.summary.files(summary.files)
  
  lambda <- validate.lambda.summaryData(summary.files, lambda)
  
  sf <- load.summary.files(summary.files, lambda, sel.snps)
  
  ref.allele <- extract.reference.allele(sf$stat)
  
  conf.snps <- extract.conflictive.snps(sf$stat, ref.allele)
  
  rcs <- remove.conflictive.snps(sf$stat, ref.allele, conf.snps)
  
  meta.stat <- merge.stat(rcs$stat, rcs$ref.allele, conf.snps, sf$lambda)
  
  meta.stat
  
}


