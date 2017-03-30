
extract.reference.allele <- function(stat){
  
  msg <- paste("Extracting allele information:", date())
  message(msg)
  
  snps <- NULL
  nstudy <- length(stat)
  for(i in 1:nstudy){
    snps <- unique(c(snps, stat[[i]][, 'SNP']))
  }
  
  nsnp <- length(snps)
  
  RefAllele <- rep(NA, nsnp)
  EffectAllele <- rep(NA, nsnp)
  names(RefAllele) <- snps
  names(EffectAllele) <- snps
  
  # use the same strategy as in recover.stat() to determine which allele is used as the reference allele when merging data or creating scores from multiple studies
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    
    if(any(is.na(RefAllele[s]))){
      s <- s[is.na(RefAllele[s])]
      RefAllele[s] <- stat[[i]][, 'RefAllele']
    }
    
    s <- stat[[i]][, 'SNP']
    if(any(is.na(EffectAllele[s]))){
      s <- s[is.na(EffectAllele[s])]
      EffectAllele[s] <- stat[[i]][, 'EffectAllele']
    }
  }
  
  list(RefAllele = RefAllele, EffectAllele = EffectAllele)
  
}
