
validate.reference <- function(reference){
  
  # validate reference
  tmp <- (c("data.frame", "matrix", 'list') %in% class(reference))
  if(!any(tmp)){
    msg <- "reference should be either an external file name or a data.frame"
    stop(msg)
  }else{
    if("matrix" %in% class(reference)){
      reference <- as.data.frame(reference)
    }
    
    if('list' %in% class(reference)){
      if(all(c('allele.info', 'ref.geno') %in% names(reference))){ # ref.does
        if(!('data.frame' %in% class(reference$allele.info))){
          msg <- 'reference$allele.info should be a data.frame'
          stop(msg)
        }
        miss.col <- setdiff(c("Chr", "SNP", "Pos", "RefAllele", "EffectAllele"), colnames(reference$allele.info))
        if(length(miss.col)==0){
          # everything is fine, do nothing
        }else{
          msg <- paste0('Columns ', paste0(miss.col, sep=', ', collapse = ''), 'are missing in reference$allele.info', sep='')
          stop(msg)
        }
      }else{
        msg <- 'reference should contain allele.info and ref.geno if it is a list'
        stop(msg)
      }
    }
  }
  
  if(reference.type(reference) %in% c('ref.geno', 'ref.does')){
    return(NULL)
  }
  
  header <- c("bed", "bim", "fam")
  tmp <- (header %in% colnames(reference))
  if(!any(tmp)){
    msg <- paste("Columns below were not found in reference:\n", paste(header[!tmp], collapse = " "))
    stop(msg)
  }
  
  reference <- reformat.reference.path(reference)
  
  tmp <- !file.exists(reference$bed)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$bed[tmp]), collapse = "\n")
    stop(msg)
  }
  
  tmp <- !file.exists(reference$bim)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$bim[tmp]), collapse = "\n")
    stop(msg)
  }
  
  tmp <- !file.exists(reference$fam)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$fam[tmp]), collapse = "\n")
    stop(msg)
  }
  
}
