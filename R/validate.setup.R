
validate.setup <- function(setup.file){
  
  ld <- try(load(setup.file), silent = TRUE)
  if(error.try(ld)){
    msg <- paste("Cannot load", setup.file)
    stop(msg)
  }
  
  if(!("setup" %in% ls())){
    msg <- paste("Cannot find R object \"setup\" in", setup.file)
    stop(msg)
  }
  
  rm(setup)
  gc()
  
}
