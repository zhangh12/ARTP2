
load.reference.allele <- function(reference, snps.in.pathway, options){
  
  msg <- paste("Loading allele information from PLINK files:", date())
  if(options$print) message(msg)
  
  if("matrix" %in% class(reference)){
    reference <- as.data.frame(reference)
  }
  
  col.class <- rep("NULL", 6)
  col.class[c(2, 5, 6)] <- "character"
  bim.files <- reference$bim
  allele.info <- NULL
  for(file in bim.files){
    bim <- read.table(file, header = FALSE, as.is = TRUE, colClasses = col.class)
    colnames(bim) <- c("SNP", "RefAllele", "EffectAllele")
    bim <- bim[bim$SNP %in% snps.in.pathway, ]
    allele.info <- rbind(allele.info, bim)
    rm(bim)
    gc()
  }
  
  if(is.null(allele.info)){
    msg <- "No SNPs were found in bim files"
    stop(msg)
  }
  
  allele.info
  
}
