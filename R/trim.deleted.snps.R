
trim.deleted.snps <- function(deleted.snps, options){
  
  if(options$tidy && nrow(deleted.snps)){
    keep.row <- which(!(deleted.snps$reason %in% c("NO_SUM_STAT", "NO_REF")))
    deleted.snps <- deleted.snps[keep.row, , drop = FALSE]
  }
  
  deleted.snps
  
}
