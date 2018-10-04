
filter.ambig.snps <- function(sum.stat, allele.info, options){
	
  msg <- paste("Removing SNPs with ambiguous alleles and missing allele frequency:", date())
  if (options$print) message(msg)

  exc.snps <- NULL
  sum.info <- sum.stat$stat
  nstudy   <- length(sum.info)	
  ref      <- allele.info[, 'RefAllele']
  eff      <- allele.info[, 'EffectAllele']
  tmp      <- ambig.snps(ref, eff)
  if (any(tmp)) {
    ambig <- (allele.info$SNP)[tmp]  
    for (k in 1:nstudy){
      dat <- sum.info[[k]]
      vec <- as.numeric(dat[, "EAF"])
      tmp <- (dat[, "SNP"] %in% ambig) & ( (!is.finite(vec)) | (vec < 0) | (vec > 1) )
      tmp[is.na(tmp)] <- TRUE
      if (any(tmp)) exc.snps <- c(exc.snps, dat[tmp, "SNP"])
    }
    exc.snps <- unique(exc.snps)
  }

  exc.snps

}


