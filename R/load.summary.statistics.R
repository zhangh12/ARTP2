
load.summary.statistics <- function(summary.files, snps.in.pathway, options){
  
  snps.in.pathway <- unique(snps.in.pathway)
  
  msg <- paste("Loading summary statistics:", date())
  if(options$print) message(msg)
  
  header <- c("SNP", "RefAllele", "EffectAllele", "BETA") # columns that must be provided by users
  opt.header <- c("P", "SE")
  
  complete.header <- c(header, opt.header, "Direction")
  
  # summary.files is a vector of file names
  stat <- list()
  lambda <- NULL
  sample.size <- list()
  fid <- 0
  snps.in.study <- NULL
  nfiles <- length(summary.files)
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
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, colClasses = col.class)
    st <- st[st$SNP %in% snps.in.pathway, , drop = FALSE]
    if(nrow(st) == 0){
      next
    }
    
    if(!("Direction" %in% colnames(st))){
      st$Direction <- paste(rep('+', length(options$sample.size[[i]])), collapse = '', sep = '')
      msg <- paste0('Direction is not found in ', summary.files[i], '. ARTP3 assumes equal sample size of SNPs in the study. ')
      warning(msg)
    }
    
    if(any(opt.header %in% colnames(st))){
      if("P" %in% colnames(st)){
        st$SE <- NA
      }else{
        st$P <- NA
      }
    }else{
      msg <- paste0("Neither SE nor P is not provided in ", summary.files[i])
      stop(msg)
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
    
    tmp <- strsplit(st$Direction, split = '')
    len <- sapply(tmp, length)
    if(length(unique(len)) > 1){
      msg <- paste0('Invalid Direction in ', summary.files[i])
      stop(msg)
    }
    
    len <- len[1]
    
    if(len != length(options$sample.size[[i]])){
      msg <- paste0('Invalid sample.size for study in ', summary.files[i])
      stop(msg)
    }
    
    fid <- fid + 1
    st$Direction <- sapply(tmp, function(x){paste(ifelse(x=='?',0,1),sep='',collapse = '')})
    #st$N <- sapply(tmp, function(x, ss = options$sample.size[[i]]){sum(ss[x != '?'])})
    rm(tmp)
    gc()
    
    stat[[fid]] <- st
    snps.in.study <- unique(c(snps.in.study, st$SNP))
    lambda[fid] <- options$lambda[i]
    sample.size[[fid]] <- options$sample.size[[i]]
    
    rm(st, len)
    gc()
  }
  
  if(length(stat) == 0){
    msg <- "No SNPs to be included in analysis"
    stop(msg)
  }
  
  snps.in.study <- unique(snps.in.study)
  list(stat = stat, snps.in.study = snps.in.study, lambda = lambda, sample.size = sample.size)
  

}

