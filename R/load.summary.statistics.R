
load.summary.statistics <- function(summary.files, snps.in.pathway, options){
  
  msg <- paste("Loading summary statistics:", date())
  if(options$print) message(msg)
  
  header <- c("SNP", "RefAllele", "EffectAllele", "BETA") # columns that must be provided by users
  opt.header <- c("P", "SE")
  
  complete.header <- c(header, opt.header, "N")
  
  # summary.files is a vector of file names
  stat <- list()
  fid <- 0
  snps.in.study <- NULL
  nfiles <- length(summary.files)
  for(i in 1:nfiles){
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, nrows = 10)
    tmp <- (header %in% colnames(st))
    if(!all(tmp)){
      msg <- paste0("Columns below were not found in ", summary.files[i], ":\n", paste(header[!tmp], collapse = " "))
      stop(msg)
    }
    
    col.class <- sapply(st, class)
    col.id <- which(colnames(st) %in% complete.header)
    col.class[-col.id] <- "NULL"
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, colClasses = col.class)
    st <- st[st$SNP %in% snps.in.pathway, , drop = FALSE]
    if(nrow(st) == 0){
      if(i < length(summary.files)){
        next
      }else{
        msg <- "No SNPs to be included in analysis"
        stop(msg)
      }
    }
    
    if(!("N" %in% colnames(st))){
      st$N <- NA
    }
    
    tmp <- (opt.header %in% colnames(st))
    if(any(tmp)){
      if("P" %in% colnames(st)){
        st$SE <- NA
      }else{
        st$P <- NA
      }
    }else{
      if(!any(tmp)){
        msg <- paste0("Neither SE nor P is not provided in ", summary.files[i])
        stop(msg)
      }
    }
    
    st <- st[, complete.header]
    
    if(nrow(st) == 0){
      next
    }
    
    dup <- duplicated(st$SNP)
    if(any(dup)){
      dup.snps <- unique(st$SNP[dup])
      msg <- paste("SNPs below are duplicated: ", paste(dup.snps, collapse = " "))
      stop(msg)
    }
    
    id.no.SE.P <- which(is.na(st$SE) & is.na(st$P))
    if(length(id.no.SE.P) > 0){
      msg <- paste("For SNPs below, neither SE nor P is not provided in", summary.files[i], ":\n", paste(st$SNP[id.no.SE.P], collapse = " "))
      stop(msg)
    }
    
    fid <- fid + 1
    stat[[fid]] <- st
    snps.in.study <- c(snps.in.study, st$SNP)
    rm(st)
    gc()
  }
  
  snps.in.study <- unique(snps.in.study)
  list(stat = stat, snps.in.study = snps.in.study)
  

}

