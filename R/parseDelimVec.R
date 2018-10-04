parseDelimVec <- function(vec, sep, ncol, numeric=0) {

  mat <- unlist(strsplit(vec, sep, fixed=TRUE))
  if (length(mat) != length(vec)*ncol) {
    stop("ERROR: check ncol or if some elements of the vector are missing delimiters")
  }
  if (numeric) mat <- as.numeric(mat)
  mat <- matrix(mat, byrow=TRUE, ncol=ncol)
  return(mat)

  mat   

} 

