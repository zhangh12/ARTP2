
artp3.chr <- function(group.setup, gene.cutpoint.setup, U, score0, V, options, cid){
  
  ngene <- length(group.setup$GeneInGroup)
  nsnp <- length(score0)
  s2 <- diag(V)
  vU <- as.vector(t(U)) # expand by row
  vV <- as.vector(V)
  rm(U, V)
  gc()
  
  method <- options$method
  nthread <- options$nthread
  nperm <- options$nperm
  ngap <- min(10000, nperm);
  nblock = nperm %/% ngap;
  seed <- options$seed + (cid - 1) * nthread * nblock
  out.dir <- options$out.dir
  id.str <- options$id.str
  file.prefix <- paste0(out.dir, "/", id.str, ".CID.", cid - 1, ".")
  
  gene.pval <- rep(1, ngene)
  arr.rank <- rep(0, length(gene.cutpoint.setup$vGeneCutPoint))
  
  METHOD <- c("adajoint_chr", "adajoint_chr", "artp3_chr")
  
  tmp <- .C(METHOD[method], as.character(file.prefix), as.integer(method), as.integer(nperm),
            as.integer(seed), as.integer(nthread), as.integer(nsnp), 
            as.integer(ngene), as.double(vU), as.double(score0), as.double(vV), 
            as.integer(group.setup$vGeneIdx), 
            as.integer(group.setup$GeneStartEnd[, "Start"]), 
            as.integer(group.setup$GeneStartEnd[, "End"]), 
            as.integer(gene.cutpoint.setup$vGeneCutPoint), 
            as.integer(gene.cutpoint.setup$GeneCutPointStartEnd[, "Start"]), 
            as.integer(gene.cutpoint.setup$GeneCutPointStartEnd[, "End"]), 
            gene.pval = as.double(gene.pval), 
            arr.rank = as.integer(arr.rank))
  
  gene.pval <- tmp$gene.pval
  names(gene.pval) <- group.setup$GeneInGroup
  
  arr.rank <- tmp$arr.rank
  
  model <- list()
  for(i in 1:length(group.setup$GeneInGroup)){
    idx <- group.setup$GeneStartEnd[i, "Start"]:group.setup$GeneStartEnd[i, "End"]
    gene.idx <- group.setup$vGeneIdx[idx]
    chi <- sort(score0[gene.idx]^2/s2[gene.idx], decreasing = TRUE)
    rs <- names(chi)
    
    cp.idx <- gene.cutpoint.setup$GeneCutPointStartEnd[i, "Start"]:gene.cutpoint.setup$GeneCutPointStartEnd[i, "End"]
    cp <- gene.cutpoint.setup$vGeneCutPoint[cp.idx]
    
    ar <- arr.rank[cp.idx]
    id <- which.min(ar)
    
    sel.snp.pvalue <- pchisq(chi[1:cp[id]], df = 1, lower.tail = FALSE)
    unadj.pvalue <- (min(ar) + 1)/(nperm + 1)
    model[[i]] <- list(sel.snp.pvalue = sel.snp.pvalue, unadj.pvalue = unadj.pvalue)
  }
  
  names(model) <- group.setup$GeneInGroup
  
  list(gene.pval = gene.pval, model = model)
  
}

