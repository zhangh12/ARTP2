
reformat.reference.path <- function(reference){
  
  if(!('list' %in% class(reference))){
    reference <- as.data.frame(reference)
  }
  
  if(reference.type(reference) %in% c('ref.geno', 'ref.does')){
    return(reference)
  }
  
  if("bed" %in% colnames(reference)){
    reference$bed <- as.character(reference$bed)
  }
  
  if("bim" %in% colnames(reference)){
    reference$bim <- as.character(reference$bim)
  }
  
  if("fam" %in% colnames(reference)){
    reference$fam <- as.character(reference$fam)
  }
  
  reference
  
}

