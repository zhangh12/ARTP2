
load.reference.geno <- function(reference, snps.in.pathway, options){
  
  msg <- paste("Loading genotypes from PLINK files:", date())
  if(options$print) message(msg)
  
  sel.subs <- options$selected.subs
  exc.subs <- options$excluded.subs
  
  if(is.null(sel.subs) && !is.null(exc.subs)){
    col.class <- rep("NULL", 6)
    col.class[2] <- "character"
    sel.subs <- read.table(reference$fam[1], header = FALSE, as.is = TRUE, colClasses = col.class)[, 1]
    sel.subs <- setdiff(sel.subs, exc.subs)
    exc.subs <- NULL
  }
  
  geno <- NULL
  for(i in 1:nrow(reference)){
    g <- read.bed(reference$bed[i], reference$bim[i], reference$fam[i], snps.in.pathway, sel.subs)
    if(!is.null(g)){
      if(is.null(geno)){
        geno <- g
      }else{
        geno <- cbind(geno, g)
      }
    }
  }
  
  geno
  
}
