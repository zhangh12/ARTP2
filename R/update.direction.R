

update.direction <- function(stat, ref.allele){
  
  RefAllele <- ref.allele$RefAllel
  EffectAllele <- ref.allele$EffectAllele
  
  snps <- names(RefAllele)
  nsnp <- length(snps)
  
  A1 <- paste(RefAllele, EffectAllele, sep = '')
  names(A1) <- snps
  
  nstudy <- length(stat)
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    A2 <- paste(stat[[i]][, 'RefAllele'], stat[[i]][, 'EffectAllele'], sep = '')
    names(A2) <- s
    id <- which(A2 != A1[s])
    if(length(id) == 0){
      next
    }
    d <- strsplit(stat[[i]][id, 'Direction'], '')
    d <- sapply(d, function(u){paste(ifelse(u %in% c('+','-'), ifelse(u == '+', '-', '+'), u), collapse = '')})
    stat[[i]][id, 'Direction'] <- d
    stat[[i]][id, 'RefAllele'] <- RefAllele[s[id]]
    stat[[i]][id, 'EffectAllele'] <- EffectAllele[s[id]]
    stat[[i]][id, 'BETA'] <- -stat[[i]][id, 'BETA']
  }
  
  stat
  
}
