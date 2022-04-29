library(bigsnpr)
obj.bigsnp <- snp_attach("tmp-data/GWAS_data_sorted_QC.rds")

G <- obj.bigsnp$genotypes

PC <- predict(readRDS("tmp-data/PCA_GWAS_data.rds"))

covar <- cbind(as.matrix(obj.bigsnp$fam[c("sex", "age")]), PC[, 1:6])

y <- obj.bigsnp$fam$CAD

ind.gwas <- which(!is.na(y) & complete.cases(covar))

gwas <- runonce::save_run(
  big_univLogReg(G, y[ind.gwas], ind.train = ind.gwas,
                       covar.train = covar[ind.gwas, ], 
                       ncores = nb_cores()),
  file = "tmp-data/GWAS_CAD.rds")

plot(gwas)

CHR <- obj.bigsnp$map$chromosome
POS <- obj.bigsnp$map$physical.pos
snp_manhattan(gwas, CHR, POS, npoints = 50e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red")

y2 <- obj.bigsnp$fam$hdl
ind.gwas2 <- which(!is.na(y2) & complete.cases(covar))
gwas2 <- big_univLinReg(G, y2[ind.gwas2], ind.train = ind.gwas2,
                        covar.train = covar[ind.gwas2, ], 
                        ncores = nb_cores())
snp_manhattan(gwas2, CHR, POS, npoints = 50e3) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red")

