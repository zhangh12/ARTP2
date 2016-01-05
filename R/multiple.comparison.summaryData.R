

multiple.comparison.summaryData <- function(summary.files, pathway, family, reference, lambda, 
                                            ncases = list(), ncontrols = list(), nsamples = list(), 
                                            options = NULL){
  
  validate.summary.input(summary.files, NULL, family, reference, lambda, ncases, ncontrols, nsamples)
  
  reference <- reformat.reference.path(reference)
  
  options <- options.setup(options, family, lambda, ncases, ncontrols, nsamples)
  
  pathway <- load.pathway.set(pathway, options)
  
  super.pathway <- create.super.pathway(pathway)
  
  setup <- summaryData.setup(summary.files, super.pathway, family, reference, lambda, 
                             ncases, ncontrols, nsamples, options)
  
  setup <- recreate.pathway(setup, pathway)
  
  generate.pathway.pvalue.stat(setup)
  
}


