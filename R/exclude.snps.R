
exclude.snps <- function(stat, excluded.regions){
  
  if(!all(c('Chr', 'Pos') %in% colnames(stat))){
    msg <- 'Columns \'Chr\' and \'Pos\' should be provided in stat'
    stop(msg)
  }
  
  if(!all(c('Chr', 'Start', 'End') %in% colnames(excluded.regions))){
    msg <- 'Columns \'Chr\', \'Start\' and \'End\' should be provided in excluded.regions'
    stop(msg)
  }
  
  nregions <- nrow(excluded.regions)
  exc.rows <- NULL
  for(j in 1:nregions){
    id <- which(stat$Chr == excluded.regions$Chr[j] & stat$Pos >= excluded.regions$Start[j] & stat$Pos <= excluded.regions$End[j])
    exc.rows <- c(exc.rows, id)
  }
  
  if(!is.null(exc.rows)){
    stat <- stat[-exc.rows, ]
  }
  
  stat
  
}

