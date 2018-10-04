
ambig.check.data <- function(data, file, vars=c("RAF", "EAF")) {

  tmp <- vars %in% colnames(data)
  if (!any(tmp)) stop(paste0("Neither RAF nor EAF is not provided in ", file)) 
   if (sum(tmp) == 1) {
      inc         <- vars[tmp]
      exc         <- vars[!tmp]
      data[, exc] <- 1 - as.numeric(makeVector(data[, inc])) 
   }
   data[, "EAF"] <- as.numeric(makeVector(data[, "EAF"]))
   data[, "RAF"] <- as.numeric(makeVector(data[, "RAF"]))

   data

}
