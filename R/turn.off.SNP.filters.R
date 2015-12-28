
turn.off.SNP.filters <- function(options){
  
  if(options$turn.off.filters){
    options$snp.miss.rate <- 1
    options$maf <- .0
    options$gene.R2 <- 1
    options$chr.R2 <- 1
    options$huge.gene.R2 <- 1
    options$huge.chr.R2 <- 1
    options$HWE.p <- 0
  }
  
  options
  
}
