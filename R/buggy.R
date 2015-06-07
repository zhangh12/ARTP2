

buggy <- function(pathway){
  setwd("/data/zhangh12/ARTP2")
  summary.files <- "DIAGRAM/data/DIAGRAM_hg19_GC_1.1.txt.gz"
  #pathway <- "4k_paths/hg19_4k_paths/LEE_BMP2_TARGETS_DN.txt.xls.gz"
  #pathway <- "tmp.txt"
  
  path <- "/data/zhangh12/meta-pathway/1000genomes/20130502/plink/super_population/EUR/"
  bed <- paste(path, list.files("/data/zhangh12/meta-pathway/1000genomes/20130502/plink/super_population/EUR/",pattern = "\\.bed$"), sep = "")
  bim <- paste(path, list.files("/data/zhangh12/meta-pathway/1000genomes/20130502/plink/super_population/EUR/",pattern = "\\.bim$"), sep = "")
  fam <- paste(path, list.files("/data/zhangh12/meta-pathway/1000genomes/20130502/plink/super_population/EUR/",pattern = "\\.fam$"), sep = "")
  
  #reference <- data.frame(bed, bim, fam, stringsAsFactors = FALSE)[16, , drop = FALSE]
  reference <- data.frame(bed, bim, fam, stringsAsFactors = FALSE)
  options <- list(gene.R2 = .9, chr.R2 = .9, maf = .0099, 
                  inspect.snp.n = 2, rm.gene.subset = FALSE, 
                  print = TRUE)
  res <- pathway.summaryData(summary.files, pathway, reference, options)
  res
}

