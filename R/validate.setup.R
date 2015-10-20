
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
  
  setup$options$only.setup <- NULL
  setup$options$save.setup <- NULL
  setup$options$path.setup <- NULL
  
  setup
  
}
