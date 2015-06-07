
pathway.summaryData <- function(summary.files, pathway, reference, options = NULL){
  
  setup <- summaryData.setup(summary.files, pathway, reference, options)
  
  if(options$only.setup){
    return(setup)
  }
  
  test <- summaryData.test(setup)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       options = setup$options, setup.timing = setup$timing)
  
}



