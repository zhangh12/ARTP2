
validate.summary.input <- function(summary.files, pathway, family, reference, lambda, 
                                   ncases, ncontrols, nsamples){
  
  validate.summary.files(summary.files)
  
  validate.pathway.definition(pathway)
  
  validate.family(family)
  
  validate.reference(reference)
  
  validate.lambda.summaryData(summary.files, lambda)
  
  validate.sample.size(family, lambda, ncases, ncontrols, nsamples)
  
}

