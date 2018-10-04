expand_pathway <- function(pathway, reference) {

  tmp   <- grepl("-", pathway[, "SNP"], fixed=TRUE)
  if (!any(tmp)) return(pathway)

  part1 <- pathway[tmp, , drop=FALSE]
  part2 <- pathway[!tmp, , drop=FALSE]   
  if (!nrow(part2)) part2 <- NULL

  # For part1, the only ones we want are of the form a-b, where a and b are integers.
  # There could be snp names such as 1-124234-I, so be careful
  mat <- matrix(data=NA, nrow=nrow(part1), ncol=2)
  vec <- makeVector(part1[, "SNP"])
  for (i in 1:nrow(part1)) {
    tmp <- as.numeric(strsplit(vec[i], "-")[[1]])
    if ((length(tmp) == 2) && all(is.finite(tmp))) mat[i, ] <- tmp 
  }
  tmp <- is.na(mat[, 1])
  if (any(tmp)) {
    part2 <- rbind(part2, part1[tmp, , drop=FALSE])
    part1 <- part1[!tmp, , drop=FALSE]
    if (!nrow(part1)) return(pathway) 
    mat   <- mat[!tmp, , drop=FALSE]
  }

  # Reference must have bim files
  bimf <- try(makeVector(reference[, "bim"]), silent=TRUE)
  if ("try-error" %in% class(bimf)) stop("reference must be a data frame with columns bed, bim, fam")

  chrs <- makeVector(part1[, "Chr"])
  ret  <- part2
 
  # Loop over each bim file, remember that the bim files can contain different chromosomes
  # We only need the chr, snp name and location columns
  col.class          <- rep("NULL", 6)
  col.class[2]       <- "character"
  col.class[c(1, 4)] <- 'integer'
  nc1                <- ncol(part1)
  cols1              <- colnames(part1)
  for (f in bimf) {
    bim  <- try(read.table(f, header=0, colClasses=col.class, as.is=TRUE), silent=TRUE)
    if ("try-error" %in% class(bim)) {
      msg <- paste0('Cannot load ', f)
      stop(msg)
    }
    bim.chrs <- unique(bim[, 1])
    tmp      <- chrs %in% bim.chrs 
    chrs2    <- chrs[tmp]
    m        <- length(chrs2)
    if (!m) next
    mat2     <- mat[tmp, , drop=FALSE]
    part1.2  <- part1[tmp, , drop=FALSE]
    locs     <- as.numeric(bim[, 3])
    # For each chr, get all snps in the given range
    for (j in 1:m) {
      tmp1 <- (bim[, 1] %in% chrs2[j]) & (locs >= mat2[j, 1]) & (locs <= mat2[j, 2])
      tmp1[is.na(tmp1)] <- FALSE
      m2 <- sum(tmp1)
      if (m2) {
        tmp <- matrix(data=makeVector(part1.2[j, ]), nrow=m2, ncol=nc1, byrow=TRUE) 
        colnames(tmp) <- cols1
        tmp[, "SNP"]  <- bim[tmp1, 2]
        ret           <- unique(rbind(ret, tmp))
      }
    }
  }
  
  ret <- as.data.frame(ret, stringsAsFactors=FALSE)
  for(i in 1:ncol(ret)){
    if('factor' %in% class(ret[, i])){
      ret[, i] <- as.character(ret[, i])
    }
  }
  
  ret
  
}
