
validate.summary.input <- function(summary.files, pathway, reference, lambda, sample.size){
  
  # validate summary.files
  if(!is.vector(summary.files)){
    msg <- "Invalid summary.files"
    stop(msg)
    
    tmp <- !file.exists(summary.files)
    if(any(tmp)){
      msg <- paste(c("File(s) below were not found: ", summary.files[tmp]), collapse = "\n")
      stop(msg)
    }
  }
  
  
  # validate pathway
  if(is.character(pathway)){
    if(length(pathway) > 1){
      msg <- "pathway should be either an external file name or a data.frame"
      stop(msg)
    }
    
    if(!file.exists(pathway)){
      msg <- paste0("File below was not found: \n", pathway)
      stop(msg)
    }
  }else{
    tmp <- (c("data.frame", "matrix") %in% class(pathway))
    if(!any(tmp)){
      msg <- "pathway should be either an external file name or a data.frame"
      stop(msg)
    }
  }
  
  # validate reference
  tmp <- (c("data.frame", "matrix") %in% class(reference))
  if(!any(tmp)){
    msg <- "pathway should be either an external file name or a data.frame"
    stop(msg)
  }else{
    if("matrix" %in% class(reference)){
      reference <- as.data.frame(reference)
    }
  }
  
  header <- c("bed", "bim", "fam")
  tmp <- (header %in% colnames(reference))
  if(!any(tmp)){
    msg <- paste("Columns below were not found in reference:\n", paste(header[!tmp], collapse = " "))
    stop(msg)
  }
  
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
  
  if(!is.vector(lambda)){
    msg <- 'lambda should be a numeric vector'
    stop(msg)
  }
  # each summary file should has one inflation factor
  if(length(summary.files) != length(lambda)){
    msg <- 'Each summary file should has one inflation factor'
    stop(msg)
  }
  
  if(!is.list(sample.size)){
    msg <- 'sample.size should be a list'
    stop(msg)
  }
  
  if(length(sample.size) != length(lambda)){
    msg <- 'Length of sample.size and lambda should be equal'
    stop(msg)
  }
  
}

