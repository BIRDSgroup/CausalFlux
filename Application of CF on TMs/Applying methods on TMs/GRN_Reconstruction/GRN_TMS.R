str_param <- function(GE_bin_data, wl){
  library(bnlearn)
  for(i in 1:ncol(GE_bin_data))
  {
    GE_bin_data[,i] <- as.factor(GE_bin_data[,i])
  }
  toy_dat_struct <- bnlearn::hc(GE_bin_data,whitelist = wl, max.iter = 100, score = "bde")
  toy_dat_param <- bn.fit(toy_dat_struct, GE_bin_data, method = "bayes")
  GRN_recons <- list(toy_dat_struct,toy_dat_param)
  names(GRN_recons) <- c("Structure learning", "Parameter learning")
  base::return(GRN_recons)}




#########################################################################
############################################ For TM1
#########################################################################
wl <- data.frame(from = c("A","B","B","X","E","Y","Z"), to = c("C","C","D","B","D","A","A"))

#### If exchange rate is 3.2
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/new_GE_data/grn_3.2/")
gge <- read.csv("Bin_GE_TM1_3.2.csv", header = TRUE)
#### If exchange rate is 320
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/new_GE_data/grn_320/")
gge <- read.csv("Bin_GE_TM1_320.csv", header = TRUE)
#### If exchange rate is 3200
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/new_GE_data/grn_3200/")
gge <- read.csv("Bin_GE_TM1_3200.csv", header = TRUE)

cc <- str_param(gge, wl)


setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/new_GE_data/grn_3200/")
saveRDS(cc$`Structure learning`,"Structure_learning.rds")
saveRDS(cc$`Parameter learning`,"Parameter_learning.rds")

#########################################################################
############################################ For TM2
#########################################################################
wl <- data.frame(from = c("A","A","A","E","E","E","I","I","X","X"), to = c("B","C","D","F","G","H","J","K","L","M"))

#### If exchange rate is 3.2
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/new_GE_data/grn_3.2/")
gge <- read.csv("Bin_GE_TM3_3_2.csv", header = TRUE)
#### If exchange rate is 320
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/new_GE_data/grn_320/")
gge <- read.csv("Bin_GE_TM3_320.csv", header = TRUE)
#### If exchange rate is 3200
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/new_GE_data/grn_3200/")
gge <- read.csv("Bin_GE_TM3_3200.csv", header = TRUE)

cc <- str_param(gge, wl)


setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/new_GE_data/grn_3200/")
saveRDS(cc$`Structure learning`,"Structure_learning.rds")
saveRDS(cc$`Parameter learning`,"Parameter_learning.rds")

#########################################################################
############################################ For TM3
#########################################################################
wl <- data.frame(from = c("A","A","A","E","E","E","E","I","I","I","I","X","X","X","X"), to = c("B","C","D","C","D","F","G","G","H","J","K","J","K","L","M"))

#### If exchange rate is 3.2
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/new_GE_data/grn_3.2/")
gge <- read.csv("Bin_GE_TM4_3.2.csv", header = TRUE)
#### If exchange rate is 320
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/new_GE_data/grn_320/")
gge <- read.csv("Bin_GE_TM4_320.csv", header = TRUE)
#### If exchange rate is 3200
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/new_GE_data/grn_3200/")
gge <- read.csv("Bin_GE_TM4_3200.csv", header = TRUE)

cc <- str_param(gge, wl)

setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/new_GE_data/grn_3200/")
saveRDS(cc$`Structure learning`,"Structure_learning.rds")
saveRDS(cc$`Parameter learning`,"Parameter_learning.rds")





