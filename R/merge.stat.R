

merge.stat <- function(stat, ref.allele, conf.snps, lambda, only.meta){
  
  msg <- paste("Merging summary statistics:", date())
  message(msg)
  
  RefAllele <- ref.allele$RefAllel
  EffectAllele <- ref.allele$EffectAllele
  
  snps <- names(RefAllele)
  nsnp <- length(snps)
  
  BETA <- rep(0, nsnp)
  SE <- rep(0, nsnp)
  names(BETA) <- snps
  names(SE) <- snps
  
  nstudy <- length(stat)
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    stat[[i]]$sgn <- ifelse(stat[[i]][, 'RefAllele'] == RefAllele[s] & stat[[i]][, 'EffectAllele'] == EffectAllele[s], 1, -1)
    stat[[i]][, 'SE'] <- stat[[i]][, 'SE'] * sqrt(lambda[i])
    stat[[i]][, 'P'] <- pchisq(stat[[i]][, 'BETA']^2/stat[[i]][, 'SE']^2, df = 1, lower.tail = FALSE)
    BETA[s] <- BETA[s] + stat[[i]][, 'sgn'] * stat[[i]][, 'BETA']/stat[[i]][, 'SE']^2
    SE[s] <- SE[s] + 1/stat[[i]][, 'SE']^2
  }
  
  SE <- sqrt(1/SE)
  BETA <- BETA * SE^2
  P <- pchisq(BETA^2/SE^2, df = 1, lower.tail = FALSE)
  
  SNP <- names(BETA)
  meta.stat <- data.frame(SNP = SNP, RefAllele = RefAllele[SNP], EffectAllele = EffectAllele[SNP], BETA = BETA[SNP], SE = SE[SNP], P = P[SNP], stringsAsFactors = FALSE)
  rownames(meta.stat) <- NULL
  
  header <- c('SNP', 'RefAllele', 'EffectAllele', 'BETA', 'SE', 'P')
  meta.stat <- meta.stat[, header]
  
  if(!only.meta){
    for(i in 1:nstudy){
      stat[[i]] <- stat[[i]][, header]
      rownames(stat[[i]]) <- NULL
      colnames(stat[[i]]) <- c('SNP', paste(header[-1], 'Study', i, sep = '.'))
    }
    
    for(i in 1:nstudy){
      meta.stat <- merge(meta.stat, stat[[i]], by = 'SNP', all = TRUE)
    }
    
    if(!is.null(conf.snps)){
      meta.stat$Conflictive.Allele <- (meta.stat$SNP %in% conf.snps)
    }
  }
  
  meta.stat <- meta.stat[order(meta.stat$P), ]
  rownames(meta.stat) <- NULL
  
  meta.stat
  
  
}

