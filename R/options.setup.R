
options.setup <- function(options){
  
  opt.default <- options.default()
  
  new.opt <- intersect(names(opt.default), names(options))
  if(length(new.opt) > 0){
    for(opt in new.opt){
      opt.default[[opt]] <- options[[opt]]
    }
  }
  options <- opt.default
  
  if(is.null(options$path.setup)){
    options$path.setup <- paste(options$out.dir, "/setup.", options$id.str, ".rda", sep = "")
  }
  
  if(options$only.setup){
    options$save.setup <- TRUE
  }
  
  if(options$trim.huge.chr){
    if(options$huge.gene.R2 > options$gene.R2){
      options$huge.gene.R2 <- max(options$gene.R2 - .1, 0)
    }
    if(options$huge.chr.R2 > options$chr.R2){
      options$huge.chr.R2 <- max(options$chr.R2 - .1, 0)
    }
  }
  
  options.validation(options)
  
  if(!is.null(options$excluded.subs) && !is.null(options$selected.subs)){
    options$selected.subs <- setdiff(options$selected.subs, options$excluded.subs)
    options$excluded.subs <- NULL
  }
  
  tmp <- .C("check_nthread", nthread = as.integer(options$nthread))
  options$nthread <- tmp$nthread
  
  options
  
}
