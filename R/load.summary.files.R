
load.summary.files <- function(summary.files, options){
  
  msg <- paste("Loading summary files:", date())
  if(options$print) message(msg)
  
  header <- c('SNP', 'RefAllele', 'EffectAllele', 'BETA', 'N') # columns that must be provided by users
  opt.header <- c('P', 'SE')
  opt.header2 <- c('Chr', 'Pos')
  
  complete.header <- c(header, opt.header, opt.header2)
  
  nfiles <- length(summary.files)
  stat <- list()
  fid <- 0
  for(i in 1:nfiles){
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, nrows = 1e4)
    tmp <- (header %in% colnames(st))
    if(!all(tmp)){
      msg <- paste0("Columns below were not found in ", summary.files[i], ":\n", paste(header[!tmp], collapse = " "))
      stop(msg)
    }
    
    col.class <- sapply(st, class)
    col.id <- which(colnames(st) %in% complete.header)
    col.class[-col.id] <- "NULL"
    col.class[c('SNP', 'RefAllele', 'EffectAllele')] <- 'character'
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, colClasses = col.class)
    if(nrow(st) == 0){
      next
    }
    
    if(!any(opt.header %in% colnames(st))){
      msg <- paste0("Neither SE nor P is not provided in ", summary.files[i])
      stop(msg)
    }
    
    
    if(!('P' %in% colnames(st))){
      st$P <- NA
    }
    
    if(!('SE' %in% colnames(st))){
      st$SE <- NA
    }
    
    if(!('Chr' %in% colnames(st))){
      st$Chr <- NA
    }
    
    if(!('Pos' %in% colnames(st))){
      st$Pos <- NA
    }
    
    st <- st[, complete.header]
    
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
    
    st$RefAllele <- toupper(st$RefAllele)
    st$EffectAllele <- toupper(st$EffectAllele)
    
    fid <- fid + 1
    rownames(st) <- st$SNP
    stat[[fid]] <- st
    rm(st)
    gc()
    
  }
  
  if(length(stat) == 0){
    msg <- "No SNPs to be included in analysis"
    stop(msg)
  }
  
  stat
  
}

