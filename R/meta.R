
meta <- function(summary.files, lambda = rep(1.0, length(summary.files))){
  
  validate.summary.files(summary.files)
  
  validate.lambda.summaryData(summary.files, lambda)
  
  stat <- load.summary.files(summary.files, list(print = TRUE))
  
  meta.stat <- merge.stat(stat, lambda)
  
  meta.stat
  
}
