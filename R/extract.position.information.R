extract.position.information <- function(stat){
  
  msg <- paste("Extracting SNP position information:", date())
  message(msg)
  
  pos.info <- NULL
  nstudy <- length(stat)
  for(i in 1:nstudy){
    if(all(c('SNP', 'Chr', 'Pos') %in% colnames(stat[[i]]))){
      pos.info <- rbind(pos.info, stat[[i]][, c('SNP', 'Chr', 'Pos')])
    }
  }
  
  pos.info <- pos.info[!duplicated(pos.info$SNP), ]
  rownames(pos.info) <- NULL
  pos.info
  
}
