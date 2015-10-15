
options.default <- function(){
  
  opt.default <- list(out.dir = getwd(), 
                      id.str = "PID", 
                      method = 3, 
                      nperm = 1E5, 
                      snp.miss.rate = .05, 
                      maf = .05, 
                      HWE.p = 1E-5, 
                      chr.R2 = .95, 
                      gene.R2 = .95, 
                      gene.miss.rate = 1.0, 
                      huge.gene = 1000, 
                      rm.gene.subset = TRUE, 
                      delete = TRUE, 
                      print = TRUE, 
                      tidy = TRUE, 
                      save.setup = TRUE, 
                      path.setup = NULL, 
                      only.setup = FALSE, 
                      keep.geno = FALSE, 
                      seed = 1, 
                      nthread = detectCores(), 
                      excluded.snps = NULL, 
                      selected.snps = NULL, 
                      excluded.subs = NULL, 
                      selected.subs = NULL, 
                      excluded.genes = NULL, 
                      inspect.snp.n = 5, 
                      inspect.snp.percent = 0, 
                      inspect.gene.n = 10, 
                      inspect.gene.percent = .05, 
                      trim.huge.chr = TRUE, 
                      huge.chr.size = 2000, 
                      huge.gene.R2 = .85, 
                      huge.chr.R2 = .85)
  
  opt.default
  
}
