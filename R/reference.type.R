
reference.type <- function(reference){
  
  if(('list' %in% class(reference)) & all(c('allele.info', 'ref.geno') %in% names(reference))){
    return('ref.does')
  }
  
  header <- c("bed", "bim", "fam")
  tmp <- (header %in% colnames(reference))
  if(is.data.frame(reference) && !any(tmp)){
    return('ref.geno')
  }else{
    return('ref.files')
  }
  
}

