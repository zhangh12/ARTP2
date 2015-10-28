
SNPs.files <- function(geno.files, pathway){
  
  SNP <- NULL
  file <- NULL
  col <- NULL
  ncol <- NULL
  for(f in geno.files){
    rs <- scan(f, what='character',nlines = 1, quiet = TRUE)
    nrs <- length(rs)
    SNP <- c(SNP, rs)
    file <- c(file, rep(f, nrs))
    col <- c(col, 1:length(rs))
    ncol <- c(ncol, rep(nrs, nrs))
  }
  
  sf <- data.frame(SNP, file, col, ncol, stringsAsFactors = FALSE)
  sf <- sf[sf$SNP %in% pathway$SNP, ]
  sf <- sf[!duplicated(sf$SNP), ]
  rownames(sf) <- NULL
  
  sf
  
}

