library(bigsnpr)
library(bigstatsr)

plink <- ("tmp-data/plink.exe")

system(glue::glue(
  "{plink} --bfile tmp-data/GWAS_data",
  " --make-bed --out tmp-data/GWAS_data_sorted"
))

bedfile2 <- snp_plinkQC(plink, "tmp-data/GWAS_data_sorted")

(rds <- snp_readBed2(bedfile2, ncores = nb_cores()))

obj.bigsnp <- snp_attach(rds)
str(obj.bigsnp, max.level = 2)

clinical <- bigreadr::fread2("tmp-data/GWAS_clinical.csv")

pheno <- clinical[match(obj.bigsnp$fam$family.ID, clinical$FamID), ]

stopifnot(all.equal(obj.bigsnp$fam$sex, pheno$sex))

obj.bigsnp$fam <- cbind(obj.bigsnp$fam, pheno[-c(1, 3)])

G <- obj.bigsnp$genotypes
counts <- big_counts(G)
counts[, 1:8]

hist(nbNA <- counts[4, ])

G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = nb_cores())
big_counts(G2, ind.col = 1:8)

big_counts(G, ind.col = 1:8)

G$code256

G2$code256

obj.bigsnp$genotypes <- G2

snp_save(obj.bigsnp)














