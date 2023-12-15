# load libraries
library(rio)

# import data
data <- import("ukb38321.txt")
data2 <- import("ukb40265.dta")
data2 <- data2[match(data$n_eid, data2$n_eid), ]

# variables
FVC <- "n_20151_0_0"
FEV1 <- "n_20150_0_0"
FEV1FVC <- "n_20258_0_0"
Smoking <- "n_20116_0_0"
Height <- "n_50_0_0"
Smoking2 <- "n_20116_2_0"
meanIMT2 <- c("n_22671_2_0", "n_22674_2_0", "n_22677_2_0", "n_22680_2_0")

# calculate mean imt over 4 measurements
meanIMT2 <- data2[, meanIMT2]
meanIMT2[which(rowMeans(is.na(meanIMT2)) > 1), ] <- NA
meanIMT2 <- rowMeans(meanIMT2, na.rm = T)

# calculate current and never smokers
Smoking <- data[, Smoking]
Current <- ifelse(Smoking == "Current", "Yes", "No")
Smoking2 <- data[, Smoking2]
Current2 <- ifelse(Smoking2 == "Current", "Yes", "No")

# create phenofile
names <- c(FVC, FEV1, FEV1FVC, Height)
pheno <- data[, c("n_eid", "n_eid", names)]
colnames(pheno) <- c("#FID", "IID", "FVC", "FEV1", "FEV1FVC", "Height")
pheno$Current <- Current
pheno$meanIMT2 <- meanIMT2
export(pheno, "phenofile.tsv", na = "NA", quote = F)

# create neversmokers phenofile
pheno_never <- pheno[which(Smoking == "Never" & Smoking2 %in% c("", "Never")), ]
export(pheno_never, "phenofile_neversmokers.tsv", na = "NA", quote = F)
