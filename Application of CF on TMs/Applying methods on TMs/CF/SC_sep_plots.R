#############################################################################################################
#############################################################################################################
#### Functions before 

cor_dat_df <- function(df1, kv, ex){
  
  # df1$CFMTR_FBA <- df2[[1]]
  # df1$CFMTR_FVA_min <- df2[[2]]
  # df1$CFMTR_FVA_max <- df2[[3]]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
  si <- which(colnames(df1) %in% x1)
  
  
  sc_val_df <- data.frame()
  for(i in 1:length(kv)){
    sc_val <- c()
    for(j in 1:length(si)){
      s <- cor.test(df1[df1$Settings==kv[i],1], df1[df1$Settings==kv[i],si[j]], method = "spearman")
      sc_val <- c(sc_val, s$estimate) 
    }
    sc_val_df <- rbind(sc_val_df, sc_val)
  }
  colnames(sc_val_df) <- x1
  sc_val_df$Settings <- kv_vec
  sc_val_df$Exchange <- rep(ex,nrow(sc_val_df))
  
  return(sc_val_df)}





#######################################################################################################
######################   TM1 
# Per_0"

# Obtain_TM1_overall <- function(perc){
#   setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TM1/new")
#   ###############  3.2
#   act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_1.csv", header = TRUE)
#   
#   tm <- "TM1"
#   p <- perc
#   ex <- 3.2
#   
#   kv <- c("WT","B","E","Z","A","X")
#   cfmtr_tm1_3.2_per0 <- get_overall_df(tm,p,ex,kv)
#   kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
#   
#   cfmtr_tm1_3.2_per0_corr_data <- cor_dat_df(act_3.2,cfmtr_tm1_3.2_per0, kv_vec, 3.2)
#   ############### 320
#   act_320 <- read.csv("Actual_Pred_data_320_TM_1.csv", header = TRUE)
#   
#   tm <- "TM1"
#   p <- perc
#   ex <- 320
#   
#   kv <- c("WT","B","E","Z","A","X")
#   cfmtr_tm1_320_per0 <- get_overall_df(tm,p,ex,kv)
#   kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
#   
#   
#   cfmtr_tm1_320_per0_corr_data <- cor_dat_df(act_320,cfmtr_tm1_320_per0, kv_vec, 320)
#   
#   ############### 3200
#   act_3200 <- read.csv("Actual_Pred_data_3200_TM_1.csv", header = TRUE)
#   
#   tm <- "TM1"
#   p <- perc
#   ex <- 3200
#   
#   kv <- c("WT","B","E","Z","A","X")
#   cfmtr_tm1_3200_per0 <- get_overall_df(tm,p,ex,kv)
#   kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
#   cfmtr_tm1_3200_per0_corr_data <- cor_dat_df(act_3200,cfmtr_tm1_3200_per0, kv_vec, 3200)
#   
#   
#   og <- rbind(cfmtr_tm1_3.2_per0_corr_data,cfmtr_tm1_320_per0_corr_data,cfmtr_tm1_3200_per0_corr_data)
#   og$TMS <- rep("TM1",nrow(og))
#   
#   
# return(og)}


  
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1/")
  ###############  3.2
  act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_1.csv", header = TRUE)
  
  # tm <- "TM1"
  # #p <- perc
  # ex <- 3.2
  # 
  # kv <- c("WT","B","E","Z","A","X")
  # cfmtr_tm1_3.2_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
  
  a <- cor_dat_df(act_3.2, kv_vec, 3.2)
  ############### 320
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1/")
  act_320 <- read.csv("Actual_Pred_data_320_TM_1.csv", header = TRUE)
  
  # tm <- "TM1"
  # #p <- perc
  # ex <- 320
  # 
  # kv <- c("WT","B","E","Z","A","X")
  # cfmtr_tm1_320_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
  
  
  b <- cor_dat_df(act_320, kv_vec, 320)
  
  ############### 3200
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1/")
  act_3200 <- read.csv("Actual_Pred_data_3200_TM_1.csv", header = TRUE)
  
  # tm <- "TM1"
  # #p <- perc
  # ex <- 3200
  # 
  # kv <- c("WT","B","E","Z","A","X")
  # cfmtr_tm1_3200_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
  c <- cor_dat_df(act_3200, kv_vec, 3200)
  
  
  overall_SC_cfmtr_tm1_per0_corr_data <- rbind(a,b,c)
  overall_SC_cfmtr_tm1_per0_corr_data$TMS <- rep("TM1",nrow(overall_SC_cfmtr_tm1_per0_corr_data))
  
  


# overall_SC_cfmtr_tm1_per0_corr_data <- Obtain_TM1_overall(1)

#######################################################################################################
######################   TM2 
# Per_0"

  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
  ###############  3.2
  act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_2.csv", header = TRUE)
  
  # tm <- "TM2"
  # #p <- perc
  # ex <- 3.2
  # 
  # kv <- c("WT","I","X","A")
  # cfmtr_tm1_3.2_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","I KO","X KO","A KO")
  
  p <- cor_dat_df(act_3.2, kv_vec, 3.2)
  ############### 320
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
  act_320 <- read.csv("Actual_Pred_data_320_TM_2.csv", header = TRUE)
  
  # tm <- "TM2"
  # #p <- perc
  # ex <- 320
  # 
  # kv <- c("WT","I","X","A")
  # cfmtr_tm1_320_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","I KO","X KO","A KO")
  
  
  q <- cor_dat_df(act_320, kv_vec, 320)
  
  ############### 3200
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
  act_3200 <- read.csv("Actual_Pred_data_3200_TM_2.csv", header = TRUE)
  
  # tm <- "TM2"
  # #p <- perc
  # ex <- 3200
  # 
  # kv <- c("WT","I","X","A")
  # cfmtr_tm1_3200_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","I KO","X KO","A KO")
  r <- cor_dat_df(act_3200, kv_vec, 3200)
  
  
  overall_SC_cfmtr_tm2_per0_corr_data  <- rbind(p,q,r)
  overall_SC_cfmtr_tm2_per0_corr_data$TMS <- rep("TM2",nrow(overall_SC_cfmtr_tm2_per0_corr_data))
  
  


#overall_SC_cfmtr_tm2_per0_corr_data <- Obtain_TM2_overall(2)

#######################################################################################################
######################   TM3 
# Per_0"


  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
  ###############  3.2
  act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_3.csv", header = TRUE)
  
  # tm <- "TM3"
  # #p <- perc
  # ex <- 3.2
  # 
  # kv <- c("WT","A","X")
  # cfmtr_tm1_3.2_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","A KO","X KO")
  
  t <- cor_dat_df(act_3.2, kv_vec, 3.2)
  ############### 320
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
  act_320 <- read.csv("Actual_Pred_data_320_TM_3.csv", header = TRUE)
  
  # tm <- "TM3"
  # #p <- perc
  # ex <- 320
  # 
  # kv <- c("WT","A","X")
  # cfmtr_tm1_320_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","A KO","X KO")
  
  
  u <- cor_dat_df(act_320, kv_vec, 320)
  
  ############### 3200
  setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
  act_3200 <- read.csv("Actual_Pred_data_3200_TM_3.csv", header = TRUE)
  
  # tm <- "TM3"
  # #p <- perc
  # ex <- 3200
  # 
  # kv <- c("WT","A","X")
  # cfmtr_tm1_3200_per0 <- get_overall_df(tm,ex,kv)
  kv_vec <- c("WT","A KO","X KO")
  v <- cor_dat_df(act_3200, kv_vec, 3200)
  
  
  overall_SC_cfmtr_tm3_per0_corr_data <- rbind(t,u,v)
  overall_SC_cfmtr_tm3_per0_corr_data$TMS <- rep("TM3",nrow(overall_SC_cfmtr_tm3_per0_corr_data))
  



##############################################################################
################################## Calling the functions


#overall_SC_cfmtr_tm1_per0_corr_data <- Obtain_TM1_overall("Per_0")

SC_TMS_df <- rbind(overall_SC_cfmtr_tm1_per0_corr_data,overall_SC_cfmtr_tm2_per0_corr_data,overall_SC_cfmtr_tm3_per0_corr_data)


library(tidyverse)
library(ggExtra)
library(gridExtra)

library(ggplot2)

somed <- ("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/")
  
  
dir.create(paste0(somed,"/SC_sep_plots"))


setwd(paste0(somed,"/SC_sep_plots/"))

ggplot(SC_TMS_df, aes(x = CF_S_FBA_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-S (FBA) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

#ggplot(SC_TMS_df, aes(x = CF_S_FBA_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-S (FBA) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Exchange rate")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)


ggsave("CF_S_FBA_vs_CF_S_FVA_MAX.pdf")
ggsave("CF_S_FBA_vs_CF_S_FVA_MAX.jpeg")

ggplot(SC_TMS_df, aes(x = CF_S_FVA_min_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-S (FVA min) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05,alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("CF_S_FVA_MIN_vs_CF_S_FVA_MAX.pdf")
ggsave("CF_S_FVA_MIN_vs_CF_S_FVA_MAX.jpeg")


ggplot(SC_TMS_df, aes(x = GIMME_FBA_all, y = GIMME_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FBA) and Actual")+ylab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("GIMME_FBA_vs_GIMME_FVA_MAX.pdf")
ggsave("GIMME_FBA_vs_GIMME_FVA_MAX.jpeg")



ggplot(SC_TMS_df, aes(x = GIMME_FVA_min_all, y = GIMME_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA min) and Actual")+ylab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("GIMME_FVA_min_vs_GIMME_FVA_MAX.pdf")
ggsave("GIMME_FVA_min_vs_GIMME_FVA_MAX.jpeg")


ggplot(SC_TMS_df, aes(x = CF_M_FBA_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-MTR (FBA) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("CF_MTR_FBA_vs_CF_MTR_FVA_MAX.pdf")
ggsave("CF_MTR_FBA_vs_CF_MTR_FVA_MAX.jpeg")

ggplot(SC_TMS_df, aes(x = CF_M_FVA_min_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-MTR (FVA min) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("CF_MTR_FVA_MIN_vs_CF_MTR_FVA_MAX.pdf")
ggsave("CF_MTR_FVA_MIN_vs_CF_MTR_FVA_MAX.jpeg")

########################## COMPARING CF-S WITH OTHER METHODS

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)


ggsave("TRIMER_all_vs_CF_S_FVA_MAX.pdf")
ggsave("TRIMER_all_vs_CF_S_FVA_MAX.jpeg")



ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05,alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05,alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)


ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX.pdf")
ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX.jpeg")


########################## COMPARING CF-MTR WITH OTHER METHODS

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)


ggsave("TRIMER_all_vs_CF_MTR_FVA_MAX.pdf")
ggsave("TRIMER_all_vs_CF_MTR_FVA_MAX.jpeg")



ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)

ggsave("GIMME_FVA_MAX_vs_CF_MTR_FVA_MAX.pdf")
ggsave("GIMME_FVA_MAX_vs_CF_MTR_FVA_MAX.jpeg")




##########################################################

# library(GGally)
# library(ggplot2)
# 
# 
# 
# 
# #g_ <- ggpairs(SC_TMS_df ,                # Data frame
# #              columns = 2:11,        # Columns
# #              aes(color = factor(TMS),  # Color by group (cat. variable)
# #                  alpha = 0.5),lower = list(continuous=wrap("points", position=position_jitter(height=0.05, width=0.05)))) + theme(legend.position = "top")+labs(color = "TMS")
# 
# g <- ggpairs(SC_TMS_df,                 # Data frame
#              columns = 2:11,        # Columns
#              aes(color = factor(TMS),  # Color by group (cat. variable)
#                  alpha = 0.5),upper  = list(continuous = "blank"),diag  = list(continuous = "blankDiag"),lower = list(continuous=wrap("points", position=position_jitter(height=0.05, width=0.05)))) + theme(legend.position = "top")+labs(color = "TMS")
# 
# # gpairs_lower <- function(g){
# #   g$plots <- g$plots[-(1:g$nrow)]
# #   g$yAxisLabels <- g$yAxisLabels[-1]
# #   g$nrow <- g$nrow -1
# # 
# #   g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
# #   g$xAxisLabels <- g$xAxisLabels[-g$ncol]
# #   g$ncol <- g$ncol - 1
# # 
# #   g
# # }
# # 
# # gpairs_lower(g)+geom_abline()
# 
# setwd(paste0(somed,"/SC_sep_plots/"))
# 
# 
# g+geom_abline()+ggtitle("Spearman correlation (SC) comaprisons between various methods")
# ggsave("Overall_pairs.pdf")
# ggsave("Overall_pairs.jpeg")





#################################### Calculating the fraction of cases in which CF-S > TRIMER
for(i in 1:nrow(SC_TMS_df)){
  bt[i] <- SC_TMS_df$CF_M_FVA_max_all[i] >= SC_TMS_df$CF_M_FBA_all[i]
}
#bt <- SC_TMS_df$CF_S_FVA_max_all > SC_TMS_df$GIMME_FVA_max_all
bt_val <- sum(bt) / length(SC_TMS_df$CF_M_FVA_max_all)
bt_val



















