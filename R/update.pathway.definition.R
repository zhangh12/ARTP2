
update.pathway.definition <- function(pathway, exc.snps){
  
  if(length(exc.snps) == 0){
    return(pathway)
  }
  
  id <- which(!(pathway$SNP %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.pathway.definition"
    stop(msg)
  }
  pathway <- pathway[id, ]
  
  pathway
  
}
