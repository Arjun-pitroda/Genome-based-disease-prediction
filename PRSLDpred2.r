library(bigsnpr)
obj.bigsnp <- snp_attach("tmp-data/GWAS_data_sorted_QC.rds")
G <- obj.bigsnp$genotypes
NCORES <- nb_cores()
map <- dplyr::transmute(obj.bigsnp$map,
                        chr = chromosome, pos = physical.pos,
                        a0 = allele2, a1 = allele1)

txt <- "tmp-data/sumstats_CAD.txt"

writeLines(readLines(txt, n = 3))


sumstats <- bigreadr::fread2(
  txt,
  select = c("chr", "bp_hg19", "noneffect_allele",
             "effect_allele", "beta", "se_dgc"),
  col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se"))


sumstats$n_eff <- 4 / (1 / 60801 + 1 / 123504) 

(info_snp <- tibble::as_tibble(snp_match(sumstats, map)))

maf <- snp_MAF(G, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
sd_val <- sqrt(2 * maf * (1 - maf))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05

library(ggplot2)
ggplot(dplyr::slice_sample(data.frame(sd_val, sd_ss, is_bad), n = 50e3)) +
  geom_point(aes(sd_val, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr(0.9) +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2) +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "To remove?", title = "")

df_beta <- info_snp[!is_bad, ]

for (chr in 1:22) 
  
  print(chr)
  
  corr0 <- runonce::save_run({
    
    ## indices in 'sumstats'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'G'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    
    POS2 <- snp_asGeneticPos(map$chr[ind.chr2], map$pos[ind.chr2], dir = "tmp-data")
    snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, infos.pos = POS2, 
            ncores = NCORES)
    
  }, file = paste0("tmp-data/corr_chr", chr, ".rds"))
}

for (chr in 1:22) 
  print(chr)
  
  corr0 <- readRDS(paste0("tmp-data/corr_chr", chr, ".rds"))
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, "tmp-data/corr", compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3  # file size in GB

df_beta <- dplyr::filter(df_beta, chr %in% 1:22)  # TO REMOVE (for speed here)
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))

ldsc_h2_est <- ldsc[["h2"]]

# LDpred2-inf
beta_inf <- snp_ldpred2_inf(corr, df_beta, ldsc_h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.col = df_beta$`_NUM_ID_`)
AUCBoot(pred_inf, obj.bigsnp$fam$CAD)


















