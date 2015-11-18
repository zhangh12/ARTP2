
options.setup <- function(options, lambda, ncases, ncontrols){
  
  opt.default <- options.default() # valid and default options
  
  opt.default$lambda <- lambda
  opt.default$ncases <- ncases
  opt.default$ncontrols <- ncontrols
  
  if(!is.null(ncases) && !is.null(ncontrols)){
    opt.default$sample.size <- list()
    for(i in 1:length(ncases)){
      opt.default$sample.size[[i]] <- ncases[[i]] + ncontrols[[i]]
    }
  }
  
  spec.opt <- names(options)
  if(('gene.R2' %in% spec.opt) && !('chr.R2' %in% spec.opt)){
    options$chr.R2 <- options$gene.R2
  }else{
    if(!('gene.R2' %in% spec.opt) && ('chr.R2' %in% spec.opt)){
      options$gene.R2 <- options$chr.R2
    }
  }
  
  new.opt <- intersect(names(opt.default), names(options)) # valid options specified by users
  
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
    #options$save.setup <- TRUE
  }
  
  if(options$trim.huge.chr){
    if(options$huge.gene.R2 > options$gene.R2){
      options$huge.gene.R2 <- max(options$gene.R2 - .1, 0)
    }
    
    if(options$huge.chr.R2 > options$chr.R2){
      options$huge.chr.R2 <- max(options$chr.R2 - .1, 0)
    }
  }
  
  if(options$gene.R2 > 1){
    options$gene.R2 <- 1
  }
  
  if(options$chr.R2 > 1){
    options$chr.R2 <- 1
  }
  
  if(options$huge.gene.R2 > 1){
    options$huge.gene.R2 <- 1
  }
  
  if(options$huge.chr.R2 > 1){
    options$huge.chr.R2 <- 1
  }
  
  options.validation(options)
  
  if(!is.null(options$excluded.regions)){ # format of excluded.regions has been validated in options.validation
    excluded.regions <- as.data.frame(options$excluded.regions, stringsAsFactors = FALSE)
    
    header1 <- c('Chr', 'Pos', 'Radius')
    header2 <- c('Chr', 'Start', 'End')
    if(all(header1 %in% colnames(excluded.regions))){
      excluded.regions <- excluded.regions[, header1, drop = FALSE]
      excluded.regions$Start <- excluded.regions$Pos - excluded.regions$Radius
      excluded.regions$End <- excluded.regions$Pos + excluded.regions$Radius
      excluded.regions <- excluded.regions[, header2, drop = FALSE]
    }else{
      excluded.regions <- excluded.regions[, header2, drop = FALSE]
    }
    
    options$excluded.regions <- excluded.regions
    rm(excluded.regions)
  }
  
  if(!is.null(options$excluded.subs) && !is.null(options$selected.subs)){
    options$selected.subs <- setdiff(options$selected.subs, options$excluded.subs)
    options$excluded.subs <- NULL
  }
  
  tmp <- .C("check_nthread", nthread = as.integer(options$nthread))
  options$nthread <- tmp$nthread
  
  options
  
}
