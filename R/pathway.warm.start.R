

pathway.warm.start <- function(setup, nperm = NULL, lambda = 1.0){
  
  setup <- validate.setup(setup)
  
  if(!is.null(nperm)){
    setup$options$nperm <- nperm
  }
  
  tmp <- .C("check_nthread", nthread = as.integer(setup$options$nthread))
  setup$options$nthread <- tmp$nthread
  
  test <- norm.stat.test(setup, lambda)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       options = setup$options, setup.timing = setup$setup.timing)
  
}
