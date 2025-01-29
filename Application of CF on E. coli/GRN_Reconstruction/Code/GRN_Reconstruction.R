dag_wl <- read_csv("DAG_WL.csv")
GE <- readRDS("D:/work/Git_Hub_CF/Applying CF on ecoli/GRN_Reconstruction/Code/GE.rds")

library(bnlearn)
bn_TR <- hc(GE, whitelist = dag_wl,max.iter = 100, debug = T, score = "bde", maxp = 5)

setwd("/data/users/cs20d300/WORK_MAIN/Integrated_ntwk/Ecoli/TRIMER/GRN/BDE/results/TRIMER/")
saveRDS(bn_TR, "bn_TR.RDS")


bn_params_TR  <- bn.fit(bn_TR, total_ecoli_bin, method = "bayes", debug = T)
setwd("/data/users/cs20d300/WORK_MAIN/Integrated_ntwk/Ecoli/TRIMER/GRN/BDE/results/TRIMER/")
saveRDS(bn_params_TR, "bn_params_TR.RDS")









