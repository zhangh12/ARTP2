
inflation <- function(p, is.p, na.rm = FALSE){
  
  if(is.p){
    x <- qchisq(median(p, na.rm), df = 1, lower.tail = FALSE)
  }else{
    x <- median(p, na.rm)
  }
  x/qchisq(.5, df = 1)
  
}
