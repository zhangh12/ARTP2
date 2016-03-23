

create.direction <- function(ncases, ncontrols, sum.data, case.name = NULL, control.name = NULL){
  
  if(length(ncases) != length(ncontrols)){
    msg <- 'Unequal lengths of ncases and ncontrols'
    stop(msg)
  }
  
  if(!is.data.frame(sum.data)){
    msg <- 'sum.data should be a data frame'
    stop(msg)
  }
  
  if(is.null(case.name)){
    case.name <- 'N_CASES'
  }
  
  if(is.null(control.name)){
    control.name <- 'N_CONTROLS'
  }
  
  if(!(case.name %in% colnames(sum.data))){
    msg <- paste0(case.name, ' is not found in sum.data')
    stop(msg)
  }
  
  if(!(control.name %in% colnames(sum.data))){
    msg <- paste0(control.name, ' is not found in sum.data')
    stop(msg)
  }
  
  
  nstudy <- length(ncases)
  n <- 0:(2^nstudy-1)
  
  y <- t(sapply(n, function(x){as.integer(intToBits(x))[1:nstudy]}))
  
  m1 <- as.vector(y %*% ncases)
  m0 <- as.vector(y %*% ncontrols)
  
  id1 <- match(sum.data[, case.name], m1)
  id0 <- match(sum.data[, control.name], m0)
  
  idx <- which(id1 != id0)
  
  for(i in idx){
    i1 <- which(m1 == sum.data[i, case.name])
    i0 <- which(m0 == sum.data[i, control.name])
    ii <- intersect(i1, i0)
    if(length(ii) != 1){
      msg <- 'This algorithm might not work exactly on this data'
      stop(msg)
    }
    id1[i] <- ii
    id0[i] <- ii
  }
  
  dir <- apply(y, 1, function(u){paste(ifelse(u==1,'*','?'),sep='',collapse = '')})
  sum.data$Direction <- dir[id1]
  
  sum.data
  
}
