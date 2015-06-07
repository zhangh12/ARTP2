

pathway.warm.start <- function(setup.file, nperm = NULL){
  
  validate.setup(setup.file)
  
  load(setup.file)
  
  if(!is.null(nperm)){
    setup$options$nperm <- nperm
  }
  
  test <- summaryData.test(setup)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       options = setup$options, setup.timing = setup$timing)
  
}
