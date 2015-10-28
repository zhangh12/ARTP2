
validate.setup <- function(setup){
  
  if(is.character(setup)){
    ld <- try(load(setup), silent = TRUE)
    if(error.try(ld)){
      msg <- paste('Cannot load', setup)
      stop(msg)
    }
    
    if(!is.list(setup)){
      msg <- paste('Invalid components in', setup)
      stop(msg)
    }
  }
  
  
  for(i in 1:length(setup$norm.stat$V)){
    rs <- sort(names(setup$norm.stat$score0[[i]]))
    setup$norm.stat$score0[[i]] <- setup$norm.stat$score0[[i]][rs]
    setup$norm.stat$V[[i]] <- setup$norm.stat$V[[i]][rs, rs]
  }
  
  setup$options$only.setup <- NULL
  setup$options$save.setup <- NULL
  setup$options$path.setup <- NULL
  
  setup
  
}
