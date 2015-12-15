

merge.stat <- function(stat, lambda){
  
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
  
  rm.snp <- NULL
  for(s in snps){
    for(i in 1:nstudy){
      ra <- stat[[i]][s, 'RefAllele']
      ea <- stat[[i]][s, 'EffectAllele']
      if(!is.na(ra) && !is.na(ea)){
        if(is.na(RefAllele[s]) && is.na(EffectAllele[s])){
          RefAllele[s] <- ra
          EffectAllele[s] <- ea
        }else{
          if(!setequal(c(RefAllele[s], EffectAllele[s]), c(ra, ea))){
            rm.snp <- c(rm.snp, s)
          }
        }
        break
      }
    }
  }
  
  if(!is.null(rm.snp)){
    for(i in 1:nstudy){
      stat[[i]] <- stat[[i]][!(stat[[i]]$SNP %in% rm.snp), ]
    }
  }
  
  BETA <- rep(0, nsnp)
  SE <- rep(0, nsnp)
  names(BETA) <- snps
  names(SE) <- snps
  
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    stat[[i]]$sgn <- ifelse(stat[[i]][, 'RefAllele'] == RefAllele[s] & stat[[i]][, 'EffectAllele'] == EffectAllele[s], 1, -1)
    stat[[i]][, 'SE'] <- stat[[i]][, 'SE'] * sqrt(lambda[i])
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
  
  for(i in 1:nstudy){
    stat[[i]] <- stat[[i]][, header]
    stat[[i]][, 'SE'] <- stat[[i]][, 'SE'] / sqrt(lambda[i])
    rownames(stat[[i]]) <- NULL
    colnames(stat[[i]]) <- c('SNP', paste(header[-1], 'Study', i, sep = '.'))
  }
  
  for(i in 1:nstudy){
    meta.stat <- merge(meta.stat, stat[[i]], by = 'SNP', all = TRUE)
  }
  
  if(!is.null(rm.snp)){
    meta.stat$Conflictive.Allele <- (meta.stat$SNP %in% rm.snp)
  }
  
  meta.stat <- meta.stat[order(meta.stat$P), ]
  
  meta.stat
  
  
}

