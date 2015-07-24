summaryData.setup <- function(summary.files, pathway, reference, options){
  
  start.time <- date()
  
  # reformat reference to convert factor to character
  reference <- reformat.reference.path(reference)
  
  # validate the format of main inputs
  validate.summary.input(summary.files, pathway, reference)
  
  # merge and reset options
  options <- options.setup(options)
  
  # load definition of pathway
  pathway <- load.pathway.definition(pathway, options)
  
  # load and check summary statistics
  sum.stat <- load.summary.statistics(summary.files, pathway$SNP, options)
  
  # deleted snps and their reason
  deleted.snps <- data.frame(SNP = NULL, reason = NULL, comment = NULL, stringsAsFactors = FALSE)
  exc.snps <- intersect(pathway$SNP, options$excluded.snps)
  exc.snps <- setdiff(exc.snps, options$selected.snps)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "RM_BY_USER", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  
  # update with valid/available SNPs
  exc.snps <- setdiff(pathway$SNP, sum.stat$snps.in.study)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "NO_SUM_STAT", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  
  # load SNPs and their reference and effect alleles in reference genotype
  allele.info <- load.reference.allele(reference, pathway$SNP, options)
  ref.snps <- allele.info$SNP
  
  # update with valid/available SNPs
  exc.snps <- setdiff(pathway$SNP, ref.snps)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "NO_REF", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  
  # filter out SNPs with conflictive allele information in summary statistics and reference
  exc.snps <- filter.conflictive.snps(sum.stat, allele.info, options)
  
  # update with valid/available SNPs
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "CONF_ALLELE_INFO", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  
  # load genotypes in reference
  ref.geno <- load.reference.geno(reference, pathway$SNP, options)
  
  # SNP filtering based on options
  filtered.markers <- filter.reference.geno(ref.geno, pathway, options)
  
  # update with valid/available SNPs
  exc.snps <- filtered.markers$SNP
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, 
                                      reason = filtered.markers$reason, 
                                      comment = filtered.markers$comment)
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  ref.geno <- update.ref.geno(ref.geno, exc.snps)
  
  # estimate P, SE, and N if they are not provided by users
  sum.stat <- complete.sum.stat(sum.stat, ref.geno, options)
  
  # recover the summary statistics
  norm.stat <- recover.stat(sum.stat, pathway, ref.geno, allele.info, options)
  
  if(!options$keep.geno){
    rm(ref.geno)
    gc()
    ref.geno <- NULL
  }
  
  # trim the information of deleted SNPs
  deleted.snps <- trim.deleted.snps(deleted.snps, options)
  
  msg <- paste0("Setup completed: ", date())
  if(options$print) message(msg)
  
  end.time <- date()
  setup.timing <- as.integer(difftime(strptime(end.time, "%c"), strptime(start.time, "%c"), units = "secs"))
  
  setup <- list(deleted.snps = deleted.snps, options = options, 
                pathway = pathway, norm.stat = norm.stat, 
                ref.geno = ref.geno, setup.timing = setup.timing)
  
  if(options$save.setup){
    msg <- paste0("setup file has been saved at ", options$path.setup)
    message(msg)
    save(setup, file = options$path.setup)
  }
  
  setup
  
}
