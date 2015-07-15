
options.validation <- function(options){
  
  if(is.null(options$id.str)){
    warning("options$id.str is NULL")
  }
  
  if(!any(options$method %in% 1:3)){
    stop("method should be 1 (AdaJoint), 2 (AdaJoint2), or 3 (ARTP)")
  }
  
  if(!is.numeric(options$huge.gene) || options$huge.gene < 0){
    stop("huge.gene should be non-negative integer")
  }
  
  if(is.null(options$nperm)){
    stop("nperm cannot be NULL")
  }
  
  if(!is.null(options$nperm) && options$nperm < 1000){
    stop("nperm is too small")
  }
  
  if(options$snp.miss.rate > .1){
    msg <- paste0("options$snp.miss.rate = ", options$snp.miss.rate, " might be too large")
    warning(msg)
  }
  
  if(options$inspect.snp.n <= 0){
    stop("options$inspect.snp.n should be a positive integer")
  }
  
  if(options$inspect.snp.percent < 0 || options$inspect.snp.percent > 1){
    stop("option$inspect.snp.percent should be in [0, 1]")
  }
  
  if(options$inspect.gene.n <= 0){
    stop("options$inspect.gene.n should be a positive integer")
  }
  
  if(options$inspect.gene.percent < 0 || options$inspect.gene.percent > 1){
    stop("option$inspect.gene.percent should be in [0, 1]")
  }
  
  if(length(intersect(options$excluded.snps, options$selected.snps)) > 0){
    stop("Some SNPs are specified in both options$excluded.snps and options$selected.snps")
  }
  
  if(!is.null(options$excluded.subs) && !is.null(options$selected.subs)){
    if(intersect(options$excluded.subs, options$selected.subs) > 0){
      stop("Some subject IDs are specified in both options$excluded.subs and options$selected.subs")
    }
  }
  
  if(options$gene.R2 < 0 || options$gene.R2 > 1){
    stop("gene.R2 should be in [0, 1]")
  }
  
  if(options$chr.R2 < 0 || options$chr.R > 1){
    stop("chr.R2 should be in [0, 1]")
  }
  
  if(options$huge.gene.R2 < 0 || options$huge.gene.R2 > 1){
    stop("huge.gene.R2 should be in [0, 1]")
  }
  
  if(options$huge.chr.R2 < 0 || options$huge.chr.R > 1){
    stop("huge.chr.R2 should be in [0, 1]")
  }
  
}
