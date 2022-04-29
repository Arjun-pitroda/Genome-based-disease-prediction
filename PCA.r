library(bigsnpr)

plink2 <- "tmp-data/plink2.exe"

rel <- snp_plinkKINGQC(plink2, "tmp-data/GWAS_data_sorted_QC.bed",
                       thr.king = 2^-4.5, make.bed = FALSE, ncores = nb_cores())
hist(log2(rel$KINSHIP), "FD")

obj.bed <- bed("tmp-data/GWAS_data_sorted_QC.bed")
obj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, k = 20, ncores = nb_cores()),
  file = "tmp-data/PCA_GWAS_data.rds")

plot(obj.svd)

plot(obj.svd, type = "scores", scores = 1:8, coeff = 0.6)

bed.ref <- bed("tmp-data/1000G_phase3_common_norel.bed")
proj <- runonce::save_run(
  bed_projectPCA(bed.ref, obj.bed, k = 25, ncores = nb_cores()),
  file = "tmp-data/proj-to-1000G.rds")

PC.ref <- predict(proj$obj.svd.ref)
proj2 <- proj$OADP_proj
fam2 <- bigreadr::fread2(sub_bed(bed.ref$bedfile, ".fam2"))

library(ggplot2)
source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], color = fam2$`Super Population`, size = I(2)) +
    geom_point(aes(proj2[, k1], proj2[, k2]), color = "black", alpha = 0.1) +
    theme_bigstatsr(0.5) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2), color = "Ref Pop") +
    coord_equal()
}), nrow = 3)

ldist <- log(bigutilsr::dist_ogk(proj2))
lims <- bigutilsr::hist_out(ldist)$lim
hist(ldist, "FD"); abline(v = lims, col = "red") 

sum(ldist > lims[2])



