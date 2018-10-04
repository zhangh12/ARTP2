
align.ambig.study <- function(study.list, ambig.snps, ambig.eaf, ambig.ref, ambig.eff) {

  # ambig.*  inputs are for the reference
  # Each study in study.list will be aligned with the reference

  if (!length(ambig.snps)) return(study.list)
  nstudy <- length(study.list)
  names(ambig.eaf) <- ambig.snps
  names(ambig.ref) <- ambig.snps
  names(ambig.eff) <- ambig.snps

  for (k in 1:nstudy){
      dat <- study.list[[k]]
      tmp <- dat[, "SNP"] %in% ambig.snps
      if (any(tmp)) {
        dat.eaf <- makeVector(as.numeric(dat[tmp, "EAF"]))
        dat.snp <- makeVector(dat[tmp, "SNP"]) 
        eaf2    <- ambig.eaf[dat.snp]
        ref2    <- ambig.ref[dat.snp]
        eff2    <- ambig.eff[dat.snp]

        # Lets check whether the alleles match now, because we change the alleles in the next step
        a1   <- dat[tmp, "RefAllele"]
        a2   <- dat[tmp, "EffectAllele"]
        tmp2 <- ((a1 == ref2) & (a2 == eff2)) | ((a1 == eff2) & (a2 == ref2))
        tmp2[is.na(tmp2)] <- FALSE
        if (!any(tmp2)) next 
        if (!all(tmp2)) {
          dat.snp <- dat.snp[tmp2]
          tmp     <- dat[, "SNP"] %in% dat.snp
          dat.eaf <- makeVector(as.numeric(dat[tmp, "EAF"]))
          eaf2    <- ambig.eaf[dat.snp]
          ref2    <- ambig.ref[dat.snp]
          eff2    <- ambig.eff[dat.snp] 
        }
        
        # Set the alleles to match the reference set, so that another flipping is not done later
        dat[tmp, "RefAllele"]    <- ref2
        dat[tmp, "EffectAllele"] <- eff2

        tmp      <- ((dat.eaf <= 0.5) & (eaf2 > 0.5)) | ((eaf2 <= 0.5) & (dat.eaf > 0.5)) 
        err      <- is.na(tmp)
        tmp[err] <- FALSE

        if (any(tmp)) {
          # Flip these SNPs
          snps             <- dat.snp[tmp]
          tmp              <- dat[, "SNP"] %in% snps
          dat[tmp, "BETA"] <- -as.numeric(dat[tmp, "BETA"])
          dat[tmp, "EAF"]  <- dat[tmp, "RAF"]
          dat[tmp, "RAF"]  <- 1 - as.numeric(dat[tmp, "EAF"])
          study.list[[k]]  <- dat
        } 

        if (any(err)) {
          # For this error condition, set the info to missing 
          dat[err, c("BETA", "SE", "P")] <- NA
        }
      }
  }

  study.list

}
