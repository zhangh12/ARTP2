align.ambig.meta <- function(study.list) {

  # Get all the ambiguous snps and use the first study with non-missing info containing that
  #    snp as the reference
  
  nstudy <- length(study.list)
  
  # Get all ambiguous snps
  ambig <- NULL
  for (i in 1:nstudy) {
    dat <- study.list[[i]]
    ref <- makeVector(dat[, 'RefAllele'])
    eff <- makeVector(dat[, 'EffectAllele'])
    tmp <- ambig.snps(ref, eff)
    if (any(tmp)) ambig <- unique(c(ambig, dat[tmp, "SNP"]))
  }

  nambig <- length(ambig)
  if (!nambig) return(study.list)

  # Get the reference info for these ambig snps
  ambig.eaf <- rep(NA, nambig)
  ambig.ref <- rep("", nambig)
  ambig.eff <- rep("", nambig)
  names(ambig.eaf) <- ambig
  names(ambig.ref) <- ambig
  names(ambig.eff) <- ambig
  snps             <- ambig
  for (i in 1:nstudy) {
    dat <- study.list[[i]]
    eaf <- as.numeric(makeVector(dat[, 'EAF']))
    tmp <- (dat[, "SNP"] %in% snps) & (is.finite(eaf)) 
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      snps2 <- dat[tmp, "SNP"]
      ambig.eaf[snps2] <- eaf[tmp]
      ambig.ref[snps2] <- dat[tmp, "RefAllele"]
      ambig.eff[snps2] <- dat[tmp, "EffectAllele"]
      # Remove these snps from the list
      snps <- snps[!(snps %in% snps2)]
    } 
    if (!length(snps)) break
  }

  # Loop over each study and align the ambig snps with the reference
  study.list <- align.ambig.study(study.list, ambig, ambig.eaf, ambig.ref, ambig.eff)

  study.list

}