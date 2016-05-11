
gen.stat.miss <- function(resp.var, null, family, G, X, lambda, options){
  
  if(options$impute & any(is.na(G))){
    msg <- paste('Imputing missing genotype with means:', date())
    if(options$print) message(msg)
    
    for(j in 1:ncol(G)){
      id <- which(is.na(G[, j]))
      if(length(id) > 0){
        G[id, j] <- mean(G[, j], na.rm = TRUE)
      }
    }
  }
  
  formula <- paste(resp.var, "~ . -1")
  mdl0 <- glm(formula, data = null, family = family)
  ng <- ncol(G)
  nx <- ncol(X)
  
  if(family == 'binomial'){
    y0 <- mdl0$fitted.values
    A0 <- y0 * (1 - y0)
    r0 <- null[, resp.var] - y0
    tmp <- try(V <- t(G) %*% (A0 * G) - t(G) %*% (A0 * X) %*% solve(t(X) %*% (A0 * X)) %*% t(X) %*% (A0 * G), silent = TRUE)
    if(error.try(tmp)){
      msg <- "Potential existence of multicollinearity detected and ARTP2 cannot automatically deal with it right now. Please check your covariates specified in formula"
      stop(msg)
    }
    
    score0 <- as.vector(t(G) %*% r0)
  }else{
    r0 <- mdl0$residuals
    s02 <- sum(r0^2)/(length(r0) - nx)
    tmp <- try(V <- (t(G) %*% G - t(G) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% G) / s02, silent = TRUE)
    if(error.try(tmp)){
      msg <- "Potential existence of multicollinearity detected and ARTP2 cannot automatically deal with it right now. Please check your covariates specified in formula"
      stop(msg)
    }
    
    score0 <- as.vector(t(G) %*% r0 / s02)
  }
  
  if(any(is.na(G))){
    
    msg <- paste('Some genotypes are missing:', date())
    if(options$print) message(msg)
    
    V.m <- matrix(NA, ng, ng)
    score0.m <- rep(NA, ng)
    rs <- colnames(G)
    rownames(V.m) <- rs
    colnames(V.m) <- rs
    names(score0.m) <- rs
    I.ab <- matrix(NA, nx, ng)
    obs.id <- !is.na(G)
    suff.n <- t(obs.id) %*% obs.id
    
    if(family == 'binomial'){
      inv.I.aa <- solve(t(X[,,drop = FALSE]) %*% (A0 * X[,,drop = FALSE])/length(y0))
      for(k in 1:ng){
        
        id <- as.vector(which(obs.id[,k]))
        if(length(id) == 0){
          msg <- paste('All genotypes of SNP', rs[k], 'are missing')
          stop(msg)
        }
        
        if(length(id) == nrow(X)){
          mdl <- mdl0
        }else{
          mdl <- glm(formula, data = null[id,,drop=FALSE], family = 'binomial')
        }
        
        y.hat <- as.vector(mdl$fitted.values)
        res <- null[id, resp.var] - y.hat
        A <- y.hat * (1 - y.hat)
        I.ab[, k] <- t(X[id, ,drop=FALSE]) %*% (A * G[id, k]) / length(y.hat)
        V.m[k, k] <- t(G[id, k]) %*% (A * G[id, k]) / length(y.hat) - t(I.ab[, k]) %*% inv.I.aa %*% I.ab[, k]
        score0.m[k] <- sum(G[id,k] * res)
      }
        
      if(ng > 1){
        for(k in 1:(ng-1)){
          for(l in (k+1):ng){
            if(!is.na(V[k,l])){
              next
            }
            id <- as.vector(which(obs.id[, k] & obs.id[, l]))
            if(length(id) == 0){
              V.m[k, l] <- 0
              V.m[l, k] <- 0
              next
            }
            
            mdl <- glm(formula, data = null[id,,drop=FALSE], family = 'binomial')
            
            y.hat <- mdl$fitted.values
            A <- y.hat * (1 - y.hat)
            V.m[k, l] <- t(G[id, k]) %*% (A * G[id, l]) / length(y.hat) - t(I.ab[, k]) %*% inv.I.aa %*% I.ab[, l]
            V.m[l, k] <- V.m[k, l]
          }
        }
      }
      
    }else{
      inv.I.aa <- solve(t(X[,,drop=FALSE]) %*% X[,,drop=FALSE] / s02/ng)
      for(k in 1:ng){
        
        id <- as.vector(which(obs.id[,k]))
        if(length(id) == 0){
          msg <- paste('All genotypes of SNP', rs[k], 'are missing')
          stop(msg)
        }
        
        if(length(id) == nrow(X)){
          mdl <- mdl0
        }else{
          mdl <- lm(formula,data=null,subset = id)
        }
        
        res <- mdl$residuals
        s2 <- sum(res^2)/(length(res)-nx)
        I.ab[,k] <- t(X[id,,drop=FALSE]) %*% G[id,k]/s2/length(res)
        V.m[k,k] <- t(G[id,k]) %*% G[id,k] /s2/length(res) - t(I.ab[,k]) %*% inv.I.aa %*% I.ab[,k]
        score0.m[k] <- sum(G[id,k] *res)/s2
      }
      
      if(ng > 1){
        for(k in 1:(ng-1)){
          for(l in (k+1):ng){
            if(!is.na(V[k,l])){
              next
            }
            
            id <- as.vector(which(obs.id[,k] & obs.id[,l]))
            if(length(id) == 0){
              V.m[k,l] <- 0
              V.m[l,k] <- 0
              next
            }
            
            mdl <- lm(formula,data = null, subset = id)
            
            res <- mdl$residuals
            s2 <- sum(res^2)/(length(res)-nx)
            V.m[k,l] <- t(G[id,k]) %*% G[id,l] /s2/length(res) - t(I.ab[,k]) %*% inv.I.aa %*% I.ab[,l]
            V.m[l,k] <- V.m[k,l]
          }
        }
      }
    }
    
    V.m <- V.m * suff.n
    
    V[is.na(V)] <- V.m[is.na(V)]
    score0[is.na(score0)] <- score0.m[is.na(score0)]
    
  }
  
  score0 <- score0 / sqrt(nrow(X)) / sqrt(lambda)
  V <- V / nrow(X)
  names(score0) <- colnames(V)
  
  rs <- sort(names(score0))
  score0 <- score0[rs]
  V <- V[rs, rs, drop = FALSE]
  
  return(list(score0 = score0, V = V))
  
}


