
load.reference.allele <- function(reference, snps.in.pathway, options){
  
  snps.in.pathway <- unique(snps.in.pathway)
  
  msg <- paste("Loading allele information from PLINK files:", date())
  if(options$print) message(msg)
  
  if("matrix" %in% class(reference)){
    reference <- as.data.frame(reference)
  }
  
  col.class <- rep("NULL", 6)
  col.class[c(2, 5, 6)] <- "character"
  col.class[c(1, 4)] <- 'integer'
  bim.files <- reference$bim
  nfiles <- length(bim.files)
  allele.info <- NULL
  for(i in 1:nfiles){
    tmp <- try(bim <- read.table(bim.files[i], header = FALSE, as.is = TRUE, colClasses = col.class), silent = TRUE)
    if(error.try(tmp)){
      msg <- paste0('Cannot load ', bim.files[i])
      stop(msg)
    }
    
    colnames(bim) <- c("Chr", "SNP", "Pos", "RefAllele", "EffectAllele")
    bim$Reference.ID <- i
    bim <- bim[bim$SNP %in% snps.in.pathway, ]
    bim$RefAllele <- toupper(bim$RefAllele)
    bim$EffectAllele <- toupper(bim$EffectAllele)
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
