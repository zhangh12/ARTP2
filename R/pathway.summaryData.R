
pathway.summaryData <- function(summary.files, pathway, reference, lambda, ncases, ncontrols, options = NULL){
  
  setup <- summaryData.setup(summary.files, pathway, reference, lambda, ncases, ncontrols, options)
  
  if(setup$options$only.setup){
    return(setup)
  }
  
  test <- norm.stat.test(setup)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       deleted.genes = setup$deleted.genes, 
       options = setup$options, setup.timing = setup$setup.timing, 
       setup = setup)
  
}



