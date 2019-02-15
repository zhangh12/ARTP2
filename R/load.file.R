

load.file <- function(file, header, select = NULL, nrows = Inf, dec = '.'){
  
  if(!is.null(select)){
    select <- which(select != 'NULL')
  }
  
  data.table::setDF(data.table::fread(file, 
                                      header = header, 
                                      select = select, 
                                      nrows = nrows, 
                                      dec = dec, 
                                      showProgress = FALSE))
  
}
