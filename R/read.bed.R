
read.bed <- function(bed, bim, fam, sel.snps = NULL, sel.subs = NULL){
  
  col.class <- c("NULL", "character", "NULL", "NULL", "character", "character")
  bim.file <- read.table(bim, header = FALSE, as.is = TRUE, colClasses = col.class)
  nsnp <- nrow(bim.file)
  
  if(is.null(sel.snps)){
    sel.snps <- bim.file[, 1]
  }else{
    sel.snps <- unique(sel.snps)
    sel.snps <- intersect(bim.file[, 1], sel.snps)
  }
  
  sel.snp.id <- which(bim.file[, 1] %in% sel.snps)
  
  nsel <- length(sel.snp.id)
  if(nsel == 0){
    return(NULL)
  }
  
  bim.file <- bim.file[sel.snp.id, , drop = FALSE]
  sel.snps <- bim.file[, 1]
  
  col.class <- rep("NULL", 6)
  col.class[2] <- "character"
  sid <- read.table(fam, header = FALSE, as.is = TRUE, colClasses = col.class)[, 1]
  nsub <- length(sid)
  
  geno <- rep(-1, nsub * nsel)
  
  tmp <- .C("ReadBED", as.character(bed), as.integer(nsub), 
            as.integer(nsnp), as.integer(nsel), as.integer(sel.snp.id), 
            geno = as.integer(geno), PACKAGE = "ARTP2")
  
  geno <- as.data.frame(matrix(tmp$geno, nrow = nsub, byrow = FALSE))
  rownames(geno) <- sid
  colnames(geno) <- sel.snps
  
  if(!is.null(sel.subs)){
    id <- which(sel.subs %in% rownames(geno))
    if(length(id) == 0){
      msg <- paste("No subjects were left in \n", bed)
      stop(msg)
    }
    sel.subs <- sel.subs[id]
    geno <- geno[sel.subs, , drop = FALSE]
  }
  
  geno[geno == -1] <- NA
  
  geno
  
}
