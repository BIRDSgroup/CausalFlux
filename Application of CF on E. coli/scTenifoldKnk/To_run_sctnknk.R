setwd("To directory where the below mentioned files are stored")
ecoli_samples_info <- read.csv("ecoli_samples_kos_wts.csv", header = T) 

ecoli_samples_info_wt_id <- which(ecoli_samples_info$Gene.perturbated=="WT")  ## 429


load("add the directory here/raw_GE_data.RData")
View(raw_GE_data)   # 2198 4189


raw_GE_data_WT <- raw_GE_data[ecoli_samples_info_wt_id,]  # 429 4189

raw_GE_data_WT_t <- t(raw_GE_data_WT)


########################## prep work 

install.packages("scTenifoldKnk")
library(scTenifoldKnk)

x_res_1 <- scTenifoldKnk::scTenifoldKnk(countMatrix = raw_GE_data_WT_t,qc = T,gKO = "geneName", nc_nNet = 2)



