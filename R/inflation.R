
inflation <- function(p, is.p = TRUE){
  
  if(is.p){
    x <- qchisq(median(p), df = 1, lower.tail = FALSE)
  }else{
    x <- median(p)
  }
  x/qchisq(.5, df = 1)
  
}
