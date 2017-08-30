
# for SNPs named as, e.g. 1:1234, convert them to be C1P1234
# while all other pattern are kept unchanged
reformat.snps <- function(snps){
  
  snps <- as.character(snps)
  tmp1 <- suppressWarnings(which(!is.na(as.integer(gsub(':', '', snps)))))
  tmp2 <- which(substr(snps, 1, 1) != ':')
  tmp3 <- which(sapply(base::strsplit(snps, ''), tail, 1) != ':')
  tmp4 <- grep(':', snps)
  non.rs.id <- intersect(intersect(intersect(tmp1, tmp2), tmp3), tmp4)
  if(length(non.rs.id) > 0){
    snps[non.rs.id] <- paste0('C', snps[non.rs.id])
    snps[non.rs.id] <- gsub(':', 'P', snps[non.rs.id])
  }
  
  snps
  
}
