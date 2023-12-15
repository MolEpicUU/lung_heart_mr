# set options and load libraries
library(rio)
library(MendelianRandomization)
library(ieugwasr)
library(genetics.binaRies)
library(cause)
library(dplyr)

set.seed(100)

source("get_proxy.R")
source("harmonize.R")

mr.fun <- function(exposure, outcome, mv, mv.study) {
  # extract lists
  exposure.name <- exposure$name
  exposure.study <- exposure$study
  exposure.snps <- exposure$snps
  exposure.gwas <- exposure$gwas
  outcome.name <- outcome$name
  outcome.study <- outcome$study
  outcome.gwas <- outcome$gwas
  
  # import data
  exposure.snps <- import(exposure.snps)$rsid
  exposure.gwas <- import(exposure.gwas)
  outcome.gwas <- import(outcome.gwas)
  
  # format data
  exposure.gwas <- exposure.gwas[, which(colnames(exposure.gwas) != "p.value")]
  outcome.gwas <- outcome.gwas[, which(colnames(outcome.gwas) != "p.value")]
  mv <- mv[, which(colnames(mv) != "p.value")]
  
  exposure.gwas <- exposure.gwas[which(exposure.gwas$se != 0), ]
  outcome.gwas <- outcome.gwas[which(outcome.gwas$se != 0), ]
  mv <- mv[which(mv$se != 0), ]
  
  # cause
  pval_thresh <- 1e-3
  
  X <- gwas_merge(exposure.gwas, outcome.gwas, snp_name_cols = c("rsid", "rsid"),
                  beta_hat_cols = c("beta", "beta"),
                  se_cols = c("se", "se"),
                  A1_cols = c("a1", "a1"),
                  A2_cols = c("a2", "a2"))
  
  params <- est_cause_params(X, sample(X$snp, size = 1000000, replace = F))
  
  X <- X %>% mutate(pval = 2 * pnorm(abs(beta_hat_1 / seb1), lower.tail = F))
  
  X_clump <- X %>% rename(rsid = snp,
                          pval = pval) %>%
    ieugwasr::ld_clump(dat = .,
                       clump_p = pval_thresh,
                       plink_bin = genetics.binaRies::get_plink_binary(),
                       bfile = "EUR")
  keep_snps <- X_clump$rsid
  
  cause_result <- cause(X = X, variants = keep_snps, param_ests = params, pval_thresh = pval_thresh)
  
  # get proxies
  avail <- exposure.snps[which(exposure.snps %in% outcome.gwas$rsid)]
  missing <- exposure.snps[which(!exposure.snps %in% avail)]
  
  proxies <- tryCatch({
    proxies <- as.data.frame(get_ld_proxies(missing, "EUR"))
    proxies <- proxies[which(proxies$SNP_B %in% exposure.gwas$rsid & proxies$SNP_B %in% outcome.gwas$rsid), ]
    proxies <- proxies[order(abs(proxies$R), decreasing = T), ]
    proxies[which(!duplicated(proxies$SNP_A)), ]
  }, error = function(e) data.frame(SNP_A = NA, SNP_B = NA, R = NA))
  
  exposure <- exposure.gwas[match(c(avail, proxies$SNP_B), exposure.gwas$rsid), ]
  outcome <- outcome.gwas[match(c(avail, proxies$SNP_B), outcome.gwas$rsid), ]
  
  # harmonize
  colnames(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  outcome$id.exposure <- exposure.name
  outcome$id.outcome <- outcome.name
  
  harmonization <- harmonise_data(exposure, outcome)
  harmonized <- harmonization[which(harmonization$mr_keep), ]
  
  # mr
  input <- mr_input(harmonized$beta.exposure, harmonized$se.exposure, harmonized$beta.outcome, harmonized$se.outcome, snps = harmonized$SNP, exposure = exposure.name, outcome = outcome.name)
  mr <- mr_allmethods(input, method = "main")
  
  # format results
  exposure.snps.with.proxies <- exposure.snps
  exposure.snps.with.proxies[match(proxies$SNP_A, exposure.snps.with.proxies)] <- proxies$SNP_B
  harmonization <- harmonization[match(exposure.snps.with.proxies, harmonization$SNP), ]
  harmonization <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, SNP = exposure.snps, Measured = !exposure.snps %in% missing, Proxy = proxies$SNP_B[match(exposure.snps, proxies$SNP_A)], R = proxies$R[match(exposure.snps, proxies$SNP_A)], Exposure.Allele1 = harmonization$effect_allele.exposure, Exposure.Allele2 = harmonization$other_allele.exposure, Exposure.Allele1.Freq = harmonization$eaf.exposure, Outcome.Allele1 = harmonization$effect_allele.outcome, Outcome.Allele2 = harmonization$other_allele.outcome, Outcome.Allele1.Freq = harmonization$eaf.outcome, "Incompatible.Alleles" = harmonization$remove, Palindromic = harmonization$palindromic, Ambiguous = harmonization$ambiguous, Keep = harmonization$mr_keep)
  
  input <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, SNP = mr@Data@snps, Exposure.Beta = mr@Data@betaX, Exposure.SE = mr@Data@betaXse, Outcome.Beta = mr@Data@betaY, Outcome.SE = mr@Data@betaYse)
  
  result <- mr@Values
  colnames(result) <- c("Method", "Estimate", "SE", "CI.Lower", "CI.Upper", "P.value")
  result <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, result)
  
  ## mvmr
  # get proxies
  avail <- exposure.snps[which(exposure.snps %in% outcome.gwas$rsid & exposure.snps %in% mv$rsid)]
  missing <- exposure.snps[which(!exposure.snps %in% avail)]
  
  proxies <- tryCatch({
    proxies <- as.data.frame(get_ld_proxies(missing, "EUR"))
    proxies <- proxies[which(proxies$SNP_B %in% exposure.gwas$rsid & proxies$SNP_B %in% outcome.gwas$rsid & proxies$SNP_B %in% mv$rsid), ]
    proxies <- proxies[order(abs(proxies$R), decreasing = T), ]
    proxies[which(!duplicated(proxies$SNP_A)), ]
  }, error = function(e) data.frame(SNP_A = NA, SNP_B = NA, R = NA))
  
  # match
  exposure <- exposure.gwas[match(c(avail, proxies$SNP_B), exposure.gwas$rsid), ]
  outcome <- outcome.gwas[match(c(avail, proxies$SNP_B), outcome.gwas$rsid), ]
  mv <- mv[which(mv$rsid %in% c(avail, proxies$SNP_B)), ]
  
  # harmonize
  colnames(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  colnames(mv) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "id.outcome")
  
  outcome$id.outcome <- outcome.name
  outcome$id.exposure <- exposure.name
  
  mv$id.exposure <- exposure.name
  
  outcome <- rbind(outcome, mv)
  
  mv_harmonization <- harmonise_data(exposure, outcome)
  harmonized <- mv_harmonization[which(mv_harmonization$mr_keep), ]
  
  keep <- as.data.frame(table(harmonized$SNP))
  harmonized <- harmonized[which(harmonized$SNP %in% keep$Var1[which(keep$Freq == length(unique(outcome$id.outcome)))]), ]
  harmonized <- harmonized[order(harmonized$SNP), ]
  
  input1 <- harmonized[which(harmonized$id.outcome == outcome.name), ]
  input2 <- harmonized[which(!harmonized$id.outcome == outcome.name), ]
  id <- unique(input2$id.outcome)
  input2 <- lapply(id, function(x) {
    input2[which(input2$id.outcome == x), ]
  })
  input2 <- do.call(cbind, input2)
  
  # mr
  mv_input <- mr_mvinput(as.matrix(cbind(input1$beta.exposure, input2[, grepl("beta.outcome", colnames(input2))])), as.matrix(cbind(input1$se.exposure, input2[, grepl("se.outcome", colnames(input2))])), input1$beta.outcome, input1$se.outcome, exposure = c(unique(input1$id.exposure), id), snps = input1$SNP)
  mr <- mr_mvivw(mv_input)
  
  # format results
  exposure.snps.with.proxies <- exposure.snps
  exposure.snps.with.proxies[match(proxies$SNP_A, exposure.snps.with.proxies)] <- proxies$SNP_B
  id <- unique(mv_harmonization$id.outcome)
  mv_harmonization <- lapply(id, function(x) {
    mv_harmonization <- mv_harmonization[which(mv_harmonization$id.outcome == x), ]
    mv_harmonization <- mv_harmonization[match(exposure.snps.with.proxies, mv_harmonization$SNP), c("effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "remove", "palindromic", "ambiguous", "mr_keep")]
    colnames(mv_harmonization) <- c("Allele1", "Allele2", "Allele1.Freq", "Incompatible.Alleles", "Palindromic", "Ambiguous", "Keep")
    mv_harmonization
  })
  names(mv_harmonization) <- id
  mv_harmonization <- do.call(cbind, mv_harmonization)
  mv_harmonization <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, SNP = exposure.snps, Measured = !exposure.snps %in% missing, Proxy = proxies$SNP_B[match(exposure.snps, proxies$SNP_A)], R = proxies$R[match(exposure.snps, proxies$SNP_A)], Exposure.Allele1 = exposure$effect_allele.exposure[match(exposure.snps, exposure$SNP)], Exposure.Allele2 = exposure$other_allele.exposure[match(exposure.snps, exposure$SNP)], Exposure.Allele1.Freq = exposure$eaf.exposure[match(exposure.snps, exposure$SNP)], mv_harmonization)
  
  bx <- mv_input@betaX
  bxse <- mv_input@betaXse
  
  colnames(bx) <- paste0(mv_input@exposure, ".Beta")
  colnames(bxse) <-  paste0(mv_input@exposure, ".SE")
  mv_input <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, SNP = mv_input@snps, bx, bxse, Outcome.Beta = mv_input@betaY, Outcome.SE = mv_input@betaYse)
  
  mv_result <- data.frame(Exposure = exposure.name, Exposure.Study = exposure.study, Outcome = outcome.name, Outcome.Study = outcome.study, Variable = mr@Exposure, Estimate = mr@Estimate, SE = mr@StdError, CI.Lower = mr@CILower, CI.Upper = mr@CIUpper, P.value = mr@Pvalue, row.names = NULL)
  
  save(harmonization, input, result, mv_harmonization, mv_input, mv_result, cause_result, file = paste0("results/", exposure.name, "_", outcome.name, ".RData"))
}

# run lung on heart mr
smoke <- import("clean/smoke_ukbb.tsv")
height <- import("clean/height_ukbb.tsv")

smoke$id.outcome <- "Smoking"
height$id.outcome <- "Height"

id <- smoke$rsid[which(smoke$rsid %in% height$rsid)]
mv <- rbind(smoke[which(smoke$rsid %in% id), ], height[which(height$rsid %in% id), ])

exposures <- list(list(name = "FEV1", study = "UKBB", snps = "clean/fev1_snps.txt", gwas = "clean/fev1_ukbb.tsv"), 
                  list(name = "FVC", study = "UKBB", snps = "clean/fvc_snps.txt", gwas = "clean/fvc_ukbb.tsv"),
                  list(name = "FEV1FVC", study = "UKBB", snps = "clean/fev1fvc_snps.txt", gwas = "clean/fev1fvc_ukbb.tsv")
)

outcomes <- list(list(name = "IMT", study = "CHARGE", gwas = "clean/imt_charge.tsv"),
                 list(name = "Plaque", study = "CHARGE", gwas = "clean/plaque_charge.tsv"),
                 list(name = "PAD", study = "MVP", gwas = "clean/pad_mvp.tsv")
)

lapply(exposures, function(exposure) lapply(outcomes, function(outcome) mr.fun(exposure, outcome, mv, "UKBB")))

# run heart on lung mr
smoke <- import("clean/smoke_tag.tsv")
height <- import("clean/height_giant.tsv")

smoke$id.outcome <- "Smoking"
height$id.outcome <- "Height"

id <- smoke$rsid[which(smoke$rsid %in% height$rsid)]
mv <- rbind(smoke[which(smoke$rsid %in% id), ], height[which(height$rsid %in% id), ])

exposures <- list(list(name = "IMT", study = "CHARGE", snps = "clean/imt_snps.txt", gwas = "clean/imt_charge.tsv"),
                  list(name = "PAD", study = "MVP", snps = "clean/pad_snps.txt", gwas = "clean/pad_mvp.tsv")
)

outcomes <- list(list(name = "FEV1", study = "UKBB", gwas = "clean/fev1_ukbb.tsv"), 
                 list(name = "FVC", study = "UKBB", gwas = "clean/fvc_ukbb.tsv"),
                 list(name = "FEV1FVC", study = "UKBB", gwas = "clean/fev1fvc_ukbb.tsv")
)

lapply(exposures, function(exposure) lapply(outcomes, function(outcome) mr.fun(exposure, outcome, mv, c("TAG", "GIANT"))))

# run lung on heart mr in neversmokers
mv <- import("clean/height_ns_ukbb.tsv")

mv$id.outcome <- "Height"

exposures <- list(list(name = "FEV1_NS", study = "UKBB", snps = "clean/fev1_snps.txt", gwas = "clean/fev1_ns_ukbb.tsv"), 
                  list(name = "FVC_NS", study = "UKBB", snps = "clean/fvc_snps.txt", gwas = "clean/fvc_ns_ukbb.tsv"),
                  list(name = "FEV1FVC_NS", study = "UKBB", snps = "clean/fev1fvc_snps.txt", gwas = "clean/fev1fvc_ns_ukbb.tsv")
)

outcomes <- list(list(name = "meanIMT_NS", study = "UKBB", gwas = "clean/meanimt_ns_ukbb.tsv"))

lapply(exposures, function(exposure) lapply(outcomes, function(outcome) mr.fun(exposure, outcome, mv, "UKBB")))

# run heart on lung mr in neversmokers
exposures <- list(list(name = "meanIMT_NS", study = "UKBB", snps = "clean/imt_snps.txt", gwas = "clean/meanimt_ns_ukbb.tsv"))

outcomes <- list(list(name = "FEV1_NS", study = "UKBB", gwas = "clean/fev1_ns_ukbb.tsv"), 
                 list(name = "FVC_NS", study = "UKBB", gwas = "clean/fvc_ns_ukbb.tsv"),
                 list(name = "FEV1FVC_NS", study = "UKBB", gwas = "clean/fev1fvc_ns_ukbb.tsv")
)

lapply(exposures, function(exposure) lapply(outcomes, function(outcome) mr.fun(exposure, outcome, mv, "UKBB")))
