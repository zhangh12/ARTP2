
validate.outcome <- function(null, resp.var, family){
  
  resp.level <- unique(null[, resp.var])
  
  if(length(resp.level) == 2){
  	if(family == 'binomial'){
	    if(!setequal(resp.level, c(0, 1))){
	      msg <- "response variable in formula should be 0 and/or 1"
	      stop(msg)
	    }
	  }else{
	  	msg <- 'response variable in formula has only two levels, while family = \'gaussian\'. Would family be \'binomial\'?'
	  	warning(msg)
	  }
  }else{
  	if(family == 'binomial'){
  		msg <- "response variable in formula should have two levels"
  		stop(msg)
  	}
  }
  
}
