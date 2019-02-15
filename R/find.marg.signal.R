
find.marg.signal <- function(sum.stat, allele.info, options){
  
  msg <- paste("Removing SNPs close to marginal signals:", date())
  if(options$print) message(msg)
  
  stat <- sum.stat$stat
  all.snp <- sort(sum.stat$snps.in.study)
  rm(sum.stat)
  gc()
  
  nstudy <- length(stat)
  
  nsnp <- length(all.snp)
  BETA <- rep(0, nsnp)
  SE <- rep(0, nsnp)
  names(BETA) <- all.snp
  names(SE) <- all.snp
  
  RefAllele <- allele.info$RefAllele
  EffectAllele <- allele.info$EffectAllele
  names(RefAllele) <- allele.info$SNP
  names(EffectAllele) <- allele.info$SNP
  
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    stat[[i]]$sgn <- ifelse(stat[[i]][, 'RefAllele'] == RefAllele[s] & stat[[i]][, 'EffectAllele'] == EffectAllele[s], 1, -1)
    # lambda has already been adjusted for SE in complete.sum.stat()
    # SE has already been added to stat in complete.sum.stat()
    BETA[s] <- BETA[s] + stat[[i]][, 'sgn'] * stat[[i]][, 'BETA']/stat[[i]][, 'SE']^2
    SE[s] <- SE[s] + 1/stat[[i]][, 'SE']^2
  }
  
  SE <- sqrt(1/SE)
  BETA <- BETA * SE^2
  P <- pchisq(BETA^2/SE^2, df = 1, lower.tail = FALSE)
  names(P) <- names(BETA)
  
  region <- NULL
  if(any(P < options$min.marg.p)){
    idx.snp <- names(P)[P < options$min.marg.p]
    idx.snp <- allele.info[allele.info$SNP %in% idx.snp, c('Chr', 'SNP', 'Pos')]
    region <- NULL
    for(i in 1:nrow(idx.snp)){
      chr <- idx.snp$Chr[i]
      pos <- idx.snp$Pos[i]
      ai <- allele.info[allele.info$Chr == chr, c('Chr', 'SNP', 'Pos')]
      tmp <- ai[ai$Pos >= pos - options$window & ai$Pos <= pos + options$window, ]
      tmp$comment <- paste0('Close to ', idx.snp$SNP[i], ' (P = ', formatC(P[idx.snp$SNP[i]], digits=0, format='E'), ')')
      region <- rbind(region, tmp)
      rm(ai)
    }
    region <- region[!duplicated(region$SNP), ]
  }
  
  region
  
}

