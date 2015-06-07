
complete.sum.stat <- function(sum.stat, ref.geno, options){
  
  msg <- paste("Calculating P, SE, and N if not provided:", date())
  if(options$print) message(msg)
  
  nf <- length(sum.stat$stat)
  for(i in 1:nf){
    
    st <- sum.stat$stat[[i]]
    id.no.SE <- which(is.na(st$SE))
    id.no.P <- which(is.na(st$P))
    
    if(length(id.no.SE) > 0){
      z2 <- qchisq(st$P[id.no.SE], df = 1, lower.tail = FALSE)
      st$SE[id.no.SE] <- abs(st$BETA[id.no.SE]/sqrt(z2))
    }
    
    if(length(id.no.P) > 0){
      st$P[id.no.P] <- pchisq((st$BETA[id.no.P]/st$SE[id.no.P])^2, df = 1, lower.tail = FALSE)
    }
    
    if(all(is.na(st$N))){
      rg <- ref.geno[, st$SNP, drop = FALSE]
      v2 <- apply(rg, 2, var, na.rm = TRUE)
      st$N <- 1/st$SE^2/v2
    }else{
      if(any(is.na(st$N))){
        id.no.N <- which(is.na(st$N))
        id.with.N <- which(!is.na(st$N))
        rg <- ref.geno[, st$SNP[id.with.N], drop = FALSE]
        v2 <- apply(rg, 2, var, na.rm = TRUE)
        const <- median(1/st$SE[id.with.N]^2/st$N[id.with.N]/v2, na.rm = TRUE)
        rm(rg)
        gc()
        
        rg <- ref.geno[, st$SNP[id.no.N], drop = FALSE]
        v2 <- apply(rg, 2, var, na.rm = TRUE)
        st$N[id.no.N] <- 1/st$SE[id.no.N]^2/const/v2
        rm(rg)
        gc()
      }
    }
    
    sum.stat$stat[[i]] <- st
    rm(st)
    gc()
  }
  
  sum.stat
  
}

