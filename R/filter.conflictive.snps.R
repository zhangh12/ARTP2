
filter.conflictive.snps <- function(sum.stat, allele.info, options){
	
  msg <- paste("Removing SNPs with conflictive allele information:", date())
  if(options$print) message(msg)
	sum.info <- sum.stat$stat
	nstudy <- length(sum.info)
	
	exc.snps <- NULL
	for(k in 1:nstudy){
		nsnps.in.study <- nrow(sum.info[[k]])
		for(i in 1:nsnps.in.study){
			rs <- sum.info[[k]][i, "SNP"]
			nr <- which(allele.info$SNP == rs)
			ai1 <- c(sum.info[[k]][i, "RefAllele"], sum.info[[k]][i, "EffectAllele"])
			ai2 <- c(allele.info[nr, "RefAllele"], allele.info[nr, "EffectAllele"])
			if(!setequal(ai1, ai2)){
				exc.snps <- c(exc.snps, rs)
			}
		}
	}
	
	exc.snps

}
