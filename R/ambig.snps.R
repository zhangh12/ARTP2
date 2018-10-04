ambig.snps <- function(vec1, vec2) {
  # vec1    Vector of alleles
  # vec2    Other vector of alleles

  vec1 <- removeWhiteSpace(toupper(vec1))
  vec2 <- removeWhiteSpace(toupper(vec2))
  tmp  <- ((vec1 == "A") & (vec2 == "T")) |
          ((vec1 == "T") & (vec2 == "A")) | 
          ((vec1 == "C") & (vec2 == "G")) |
          ((vec1 == "G") & (vec2 == "C"))
  tmp[is.na(tmp)] <- FALSE

  tmp
} 