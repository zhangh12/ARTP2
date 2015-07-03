

pathway.warm.start <- function(setup.file, nperm = NULL, inflation.factor = NULL){
  
  validate.setup(setup.file)
  
  load(setup.file)
  
  if(!is.null(nperm)){
    setup$options$nperm <- nperm
  }
  
  if(!is.null(inflation.factor) && is.numeric(inflation.factor)){
    setup$options$inflation.factor <- inflation.factor
  }
  
  test <- norm.stat.test(setup)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       options = setup$options, setup.timing = setup$setup.timing)
  
}
