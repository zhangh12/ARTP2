
# Function to align the snps with ambiguous alleles in the summary files to the
#   reference genotype data. The aligning will be done by allele frequency.
align.ambig.snps <- function(sum.stat, allele.info, ref.geno, options){
	
  msg <- paste("Aligning SNPs with ambiguous alleles:", date())
  if (options$print) message(msg)

  sum.info <- sum.stat$stat
  nstudy   <- length(sum.info)	
  ref      <- allele.info[, 'RefAllele']
  eff      <- allele.info[, 'EffectAllele']
  tmp      <- ambig.snps(ref, eff)
  if (any(tmp)) {
    ambig      <- (allele.info$SNP)[tmp]  
    ref        <- ref[tmp]
    eff        <- eff[tmp]

    # Get the EAF for the reference snps
    ref.eaf    <- 0.5*colMeans(ref.geno[, ambig, drop=FALSE], na.rm=TRUE)

    # Loop over each study
    sum.info <- align.ambig.study(sum.info, ambig, ref.eaf, ref, eff)    
    sum.stat$stat <- sum.info
  }

  sum.stat

}


