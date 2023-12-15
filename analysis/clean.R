library(rio)

colnames <- c("rsid", "beta", "se", "p.value", "a1", "a2", "fa1")

fvc <- import("raw/run_gwas_allchr.FVC.glm.linear", format = "tsv")
fvc$rsid <- gsub("_.*$", "", fvc$ID)
fvc$a2 <- fvc$REF
fvc$a2[fvc$A1 == fvc$REF] <- fvc$ALT1[fvc$A1 == fvc$REF]
fvc$BETA[fvc$A1 == fvc$ALT] <- -fvc$BETA[fvc$A1 == fvc$ALT]
fvc <- fvc[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fvc) <- colnames
fvc$p.value[which(fvc$p.value == 1)] <- 0.99
fvc$beta <- fvc$beta / sd(fvc$beta, na.rm = T)
fvc$se <- abs(fvc$beta / qnorm(fvc$p.value / 2))
export(fvc, "clean/fvc_ukbb.tsv")

fev1 <- import("raw/run_gwas_allchr.FEV1.glm.linear", format = "tsv")
fev1$rsid <- gsub("_.*$", "", fev1$ID)
fev1$a2 <- fev1$REF
fev1$a2[fev1$A1 == fev1$REF] <- fev1$ALT1[fev1$A1 == fev1$REF]
fev1$BETA[fev1$A1 == fev1$ALT] <- -fev1$BETA[fev1$A1 == fev1$ALT]
fev1 <- fev1[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fev1) <- colnames
fev1$p.value[which(fev1$p.value == 1)] <- 0.99
fev1$beta <- fev1$beta / sd(fev1$beta, na.rm = T)
fev1$se <- abs(fev1$beta / qnorm(fev1$p.value / 2))
export(fev1, "clean/fev1_ukbb.tsv")

fev1fvc <- import("raw/run_gwas_allchr.FEV1FVC.glm.linear", format = "tsv")
fev1fvc$rsid <- gsub("_.*$", "", fev1fvc$ID)
fev1fvc$a2 <- fev1fvc$REF
fev1fvc$a2[fev1fvc$A1 == fev1fvc$REF] <- fev1fvc$ALT1[fev1fvc$A1 == fev1fvc$REF]
fev1fvc$BETA[fev1fvc$A1 == fev1fvc$ALT] <- -fev1fvc$BETA[fev1fvc$A1 == fev1fvc$ALT]
fev1fvc <- fev1fvc[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fev1fvc) <- colnames
fev1fvc$p.value[which(fev1fvc$p.value == 1)] <- 0.99
fev1fvc$beta <- fev1fvc$beta / sd(fev1fvc$beta, na.rm = T)
fev1fvc$se <- abs(fev1fvc$beta / qnorm(fev1fvc$p.value / 2))
export(fev1fvc, "clean/fev1fvc_ukbb.tsv")

pad <- import("raw/MVP.te.PAD.dbGAP.txt")
pad <- pad[, c("SNP_ID", "Effect", "SE", "PValue", "Allele1", "Allele2", "EAF")]
colnames(pad) <- colnames
pad$p.value[which(pad$p.value == 1)] <- 0.99
pad$beta <- pad$beta / sd(pad$beta, na.rm = T)
pad$se <- abs(pad$beta / qnorm(pad$p.value / 2))
export(pad, "clean/pad_mvp.tsv")

imt <- import("raw/IMT.EA.META.MAF1.HetDF4_jun.csv")
imt <- imt[, c("refid", "Effect", "StdErr", "P", "Allele1", "Allele2", "Freq1")]
colnames(imt) <- colnames
imt$p.value[which(imt$p.value == 1)] <- 0.99
imt$beta <- imt$beta / sd(imt$beta, na.rm = T)
imt$se <- abs(imt$beta / qnorm(imt$p.value / 2))
export(imt, "clean/imt_charge.tsv")

plaque <- import("raw/Plaque_meta_032218.csv")
plaque <- plaque[, c("refid", "Effect", "StdErr", "P", "Allele1", "Allele2", "Freq1")]
colnames(plaque) <- colnames
plaque$p.value[which(plaque$p.value == 1)] <- 0.99
plaque$beta <- plaque$beta / sd(plaque$beta, na.rm = T)
plaque$se <- abs(plaque$beta / qnorm(plaque$p.value / 2))
export(plaque, "clean/plaque_charge.tsv")

smoke <- import("raw/tag.evrsmk.tbl", format = "tsv")
smoke <- smoke[, c("SNP", "OR", "SE", "P", "A1", "A2", "FRQ_A")]
colnames(smoke) <- colnames
smoke$p.value[which(smoke$p.value == 1)] <- 0.99
smoke$beta <- smoke$beta / sd(smoke$beta, na.rm = T)
smoke$se <- abs(smoke$beta / qnorm(smoke$p.value / 2))
export(smoke, "clean/smoke_tag.tsv")

height <- import("raw/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt")
height <- height[, c("MarkerName", "b", "SE", "p", "Allele1", "Allele2", "Freq.Allele1.HapMapCEU")]
colnames(height) <- colnames
height$p.value[which(height$p.value == 1)] <- 0.99
height$beta <- height$beta / sd(height$beta, na.rm = T)
height$se <- abs(height$beta / qnorm(height$p.value / 2))
export(height, "clean/height_giant.tsv")

smoke <- import("raw/run_gwas_allchr.Current.glm.logistic", format = "tsv")
smoke$rsid <- gsub("_.*$", "", smoke$ID)
smoke$a2 <- smoke$REF
smoke$a2[smoke$A1 == smoke$REF] <- smoke$ALT1[smoke$A1 == smoke$REF]
smoke$BETA[smoke$A1 == smoke$ALT] <- -smoke$BETA[smoke$A1 == smoke$ALT]
smoke <- smoke[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(smoke) <- colnames
smoke$p.value[which(smoke$p.value == 1)] <- 0.99
smoke$beta <- smoke$beta / sd(smoke$beta, na.rm = T)
smoke$se <- abs(smoke$beta / qnorm(smoke$p.value / 2))
export(smoke, "clean/smoke_ukbb.tsv")

height <- import("raw/run_gwas_allchr.Height.glm.linear", format = "tsv")
height$rsid <- gsub("_.*$", "", height$ID)
height$a2 <- height$REF
height$a2[height$A1 == height$REF] <- height$ALT1[height$A1 == height$REF]
height$BETA[height$A1 == height$ALT] <- -height$BETA[height$A1 == height$ALT]
height <- height[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(height) <- colnames
height$p.value[which(height$p.value == 1)] <- 0.99
height$beta <- height$beta / sd(height$beta, na.rm = T)
height$se <- abs(height$beta / qnorm(height$p.value / 2))
export(height, "clean/height_ukbb.tsv")

fvc <- import("raw/run_gwas_allchr.FVC_NS.glm.linear", format = "tsv")
fvc$rsid <- gsub("_.*$", "", fvc$ID)
fvc$a2 <- fvc$REF
fvc$a2[fvc$A1 == fvc$REF] <- fvc$ALT1[fvc$A1 == fvc$REF]
fvc$BETA[fvc$A1 == fvc$ALT] <- -fvc$BETA[fvc$A1 == fvc$ALT]
fvc <- fvc[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fvc) <- colnames
fvc$p.value[which(fvc$p.value == 1)] <- 0.99
fvc$beta <- fvc$beta / sd(fvc$beta, na.rm = T)
fvc$se <- abs(fvc$beta / qnorm(fvc$p.value / 2))
export(fvc, "clean/fvc_ns_ukbb.tsv")

fev1 <- import("raw/run_gwas_allchr.FEV1_NS.glm.linear", format = "tsv")
fev1$rsid <- gsub("_.*$", "", fev1$ID)
fev1$a2 <- fev1$REF
fev1$a2[fev1$A1 == fev1$REF] <- fev1$ALT1[fev1$A1 == fev1$REF]
fev1$BETA[fev1$A1 == fev1$ALT] <- -fev1$BETA[fev1$A1 == fev1$ALT]
fev1 <- fev1[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fev1) <- colnames
fev1$p.value[which(fev1$p.value == 1)] <- 0.99
fev1$beta <- fev1$beta / sd(fev1$beta, na.rm = T)
fev1$se <- abs(fev1$beta / qnorm(fev1$p.value / 2))
export(fev1, "clean/fev1_ns_ukbb.tsv")

fev1fvc <- import("raw/run_gwas_allchr.FEV1FVC_NS.glm.linear", format = "tsv")
fev1fvc$rsid <- gsub("_.*$", "", fev1fvc$ID)
fev1fvc$a2 <- fev1fvc$REF
fev1fvc$a2[fev1fvc$A1 == fev1fvc$REF] <- fev1fvc$ALT1[fev1fvc$A1 == fev1fvc$REF]
fev1fvc$BETA[fev1fvc$A1 == fev1fvc$ALT] <- -fev1fvc$BETA[fev1fvc$A1 == fev1fvc$ALT]
fev1fvc <- fev1fvc[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(fev1fvc) <- colnames
fev1fvc$p.value[which(fev1fvc$p.value == 1)] <- 0.99
fev1fvc$beta <- fev1fvc$beta / sd(fev1fvc$beta, na.rm = T)
fev1fvc$se <- abs(fev1fvc$beta / qnorm(fev1fvc$p.value / 2))
export(fev1fvc, "clean/fev1fvc_ns_ukbb.tsv")

height <- import("raw/run_gwas_allchr.Height_NS.glm.linear", format = "tsv")
height$rsid <- gsub("_.*$", "", height$ID)
height$a2 <- height$REF
height$a2[height$A1 == height$REF] <- height$ALT1[height$A1 == height$REF]
height$BETA[height$A1 == height$ALT] <- -height$BETA[height$A1 == height$ALT]
height <- height[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(height) <- colnames
height$p.value[which(height$p.value == 1)] <- 0.99
height$beta <- height$beta / sd(height$beta, na.rm = T)
height$se <- abs(height$beta / qnorm(height$p.value / 2))
export(height, "clean/height_ns_ukbb.tsv")

meanimt <- import("raw/run_gwas_allchr.meanIMT_NS.glm.linear", format = "tsv")
meanimt$rsid <- gsub("_.*$", "", meanimt$ID)
meanimt$a2 <- meanimt$REF
meanimt$a2[meanimt$A1 == meanimt$REF] <- meanimt$ALT1[meanimt$A1 == meanimt$REF]
meanimt$BETA[meanimt$A1 == meanimt$ALT] <- -meanimt$BETA[meanimt$A1 == meanimt$ALT]
meanimt <- meanimt[, c("rsid", "BETA", "SE", "P", "A1", "a2", "A1_FREQ")]
colnames(meanimt) <- colnames
meanimt$p.value[which(meanimt$p.value == 1)] <- 0.99
meanimt$beta <- meanimt$beta / sd(meanimt$beta, na.rm = T)
meanimt$se <- abs(meanimt$beta / qnorm(meanimt$p.value / 2))
export(meanimt, "clean/meanimt_ns_ukbb.tsv")

fvc_snps <- import("raw/41588_2018_321_MOESM3_ESM.xlsx", which = 1, skip = 3)
fvc_snps <- fvc_snps[which(fvc_snps$Phenotype == "FVC"), ]
export(fvc_snps[, "rsid", drop = F], "clean/fvc_snps.txt")

fev1_snps <- import("raw/41588_2018_321_MOESM3_ESM.xlsx", which = 1, skip = 3)
fev1_snps <- fev1_snps[which(fev1_snps$Phenotype == "FEV1"), ]
export(fev1_snps[, "rsid", drop = F], "clean/fev1_snps.txt")

fev1fvc_snps <- import("raw/41588_2018_321_MOESM3_ESM.xlsx", which = 1, skip = 3)
fev1fvc_snps <- fev1fvc_snps[which(fev1fvc_snps$Phenotype == "FEV1/FVC"), ]
export(fev1fvc_snps[, "rsid", drop = F], "clean/fev1fvc_snps.txt")

pad_snps <- import("raw/s41591-019-0492-5.tsv")
export(pad_snps[, "rsid", drop = F], "clean/pad_snps.txt")

imt_snps <- import("raw/s41467-018-07340-5.tsv")
imt_snps <- imt_snps[c(2:9, 13:15), ]
colnames(imt_snps)[1] <- "rsid"
imt_snps$rsid <- gsub("a", "", imt_snps$rsid)
export(imt_snps[, "rsid", drop = F], "clean/imt_snps.txt")

plaque_snps <- import("raw/s41467-018-07340-5.tsv")
plaque_snps <- plaque_snps[c(11, 17:20), ]
colnames(plaque_snps)[1] <- "rsid"
plaque_snps$rsid <- gsub("b", "", plaque_snps$rsid)
export(plaque_snps[, "rsid", drop = F], "clean/plaque_snps.txt")
