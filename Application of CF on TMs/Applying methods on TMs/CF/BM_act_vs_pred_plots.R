
# obtain_df <- function(tm, p, ex, co){
#   library(readxl)
#   setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/",tm,"/",p,"/",ex,"/",co,"/"))
#   
#   y <- read.csv("FBA_to_check.csv", header = F)
#   FBA_WT <- y$V1
#   
#   z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/",tm,"/",p,"/",ex,"/",co,"/FVA_to_check.xlsx"), col_names  = F)
#   
#   z <- as.data.frame(z)
#   FVA_min_WT <- z[,2]
#   FVA_max_WT <- z[,3]
#   
#   tg_df <- data.frame(FBA_WT, FVA_min_WT, FVA_max_WT)
#   colnames(tg_df) <- c("FBA","FVA_min","FVA_max") 
#   tg_df$condition <- rep(co, nrow(tg_df))
#   
#   return(tg_df)}
# 
# get_overall_df <- function(tm, p, ex, kv){
#   ddf <- data.frame()
#   
#   for(i in 1:length(kv)){
#     wt <- obtain_df(tm,p,ex,kv[i])
#     ddf <- rbind(ddf,wt)
#   }
#   
#   
#   return(ddf)}

# obtain_df <- function(tm, ex, co){
#   library(readxl)
#   setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/",tm,"/Results/",ex,"/KO_data_",tm,"_",ex,"/",co))
#   
#   y <- read.csv("FBA_to_check.csv", header = F)
#   FBA_WT <- y$V1
#   
#   z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/",tm,"/Results/",ex,"/KO_data_",tm,"_",ex,"/",co,"/FVA_to_check.xlsx"), col_names  = F)
#   
#   z <- as.data.frame(z)
#   FVA_min_WT <- z[,2]
#   FVA_max_WT <- z[,3]
#   
#   tg_df <- data.frame(FBA_WT, FVA_min_WT, FVA_max_WT)
#   colnames(tg_df) <- c("FBA","FVA_min","FVA_max") 
#   tg_df$condition <- rep(co, nrow(tg_df))
#   
#   return(tg_df)}
# 
# 
# get_overall_df <- function(tm, ex, kv){
#   ddf <- data.frame()
#   
#   for(i in 1:length(kv)){
#     wt <- obtain_df(tm,ex,kv[i])
#     ddf <- rbind(ddf,wt)
#   }
#   
#   
#   return(ddf)}

# 
# 
# t <- 2
# 
# BM_plots <- function(t){
#   e <- c(3.2, 320, 3200)
#   
#   
#   di <- data.frame()
#   for(i in 1:length(e)){
#     setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM",t,"/"))
#     act_3.2 <- read.csv(paste0("Actual_Pred_data_",e[i],"_TM_",t,".csv"), header = TRUE)
#     
#     tm <- paste0("TM",t)
#     p <- "Per_0"
#     ex <- e[i]
#     
#     if(t == 2){
#       kv <- c("WT","I","X","A")
#     }else if(t == 1){
#       kv <- c("WT","B","E","Z","A","X")
#     }else if(t == 3){
#       kv <- c("WT","A","X")
#     }
#     
#     
#     
#     # cfmtr_tm1_3.2_per0 <- get_overall_df(tm,p,ex,kv)
#     # 
#     # act_3.2$CFMTR_FBA <- cfmtr_tm1_3.2_per0[[1]]
#     # act_3.2$CFMTR_FVA_min <- cfmtr_tm1_3.2_per0[[2]]
#     # act_3.2$CFMTR_FVA_max <- cfmtr_tm1_3.2_per0[[3]]
#     
#     d1 <- act_3.2
#     di <- rbind(di,d1)
#   }
#   
#   if(t == 1){
#     dii <- di[c(4,12,20,28,36,44,52,60,68,76,84,92,100,108,116,124,132,140),]
#     #dii <- di[c(4,12,20,28,36,44),]
#   }else if(t == 2){
#     dii <- di[c(11,22,33,44,55,66,77,88,99,110,121,132),]
#     #dii <- di[c(11,22,33,44),]
#   }else if(t == 3){
#     dii <- di[c(8,17,26,35,44,53,62,71,80),]
#   }
#   
#   
#   x <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
#   si <- which(colnames(dii) %in% x)
#   
#   
#   for(i in 1:length(si)){
#     dii[[si[i]]] <- scale(dii[[si[i]]], center = min(dii[[si[i]]]), scale = max(dii[[si[i]]]) - min(dii[[si[i]]]))
#   }
#   
#   
#   library(ggplot2)
#   library(ggpubr)
#   
#   
#   
#   somed <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/Results_1/BM_act_vs_pred/"
#   
#   dir.create(paste0(somed,"/TM",t))
#   
#   setwd(paste0(somed,"/TM",t))
#   
#   ggplot(dii, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
#   ggsave("CF_S_FVA_MAX_vs_ACTUAL.pdf")
#   ggsave("CF_S_FVA_MAX_vs_ACTUAL.jpeg")
#   
#   ggplot(dii, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
#   ggsave("TRIMER_vs_ACTUAL.pdf")
#   ggsave("TRIMER_vs_ACTUAL.jpeg")
#   
#   ggplot(dii, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
#   ggsave("GIMME_FVA_MAX_vs_ACTUAL.pdf")
#   ggsave("GIMME_FVA_MAX_vs_ACTUAL.jpeg")
#   
#   ggplot(dii, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
#   ggsave("CF_MTR_FVA_MAX_vs_ACTUAL.pdf")
#   ggsave("CF_MTR_FVA_MAX_vs_ACTUAL.jpeg")
# }

# 
# 
# BM_plots(1)
# BM_plots(2)
# BM_plots(3)
# 
# 
# 
# 


#######################################################################################################
#######################################################################################################
####### Alternate plotting  1
#######################################################################################################



#e <- c(3.2,320,3200)

all_tms_ex <- function(e){
  t <- c(1, 2, 3)
  di <- data.frame()
  
  for(i in 1:length(t)){
    setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM",t[i],"/"))
    
    #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM",t[i],"_25_75_per/"))
    
    act_3.2 <- read.csv(paste0("Actual_Pred_data_",e,"_TM_",t[i],".csv"), header = TRUE)
    
    tm <- paste0("TM",t[i])
    #p <- "Per_0"
    ex <- e
    
    if(t[i] == 2){
      kv <- c("WT","I","X","A")
    }else if(t[i] == 1){
      kv <- c("WT","B","E","Z","A","X")
    }else if(t[i] == 3){
      kv <- c("WT","A","X")
    }
    
    
    
    # cfmtr_tm1_3.2_per0 <- get_overall_df(tm,ex,kv)
    # 
    # act_3.2$CFMTR_FBA <- cfmtr_tm1_3.2_per0[[1]]
    # act_3.2$CFMTR_FVA_min <- cfmtr_tm1_3.2_per0[[2]]
    # act_3.2$CFMTR_FVA_max <- cfmtr_tm1_3.2_per0[[3]]
    
    d1 <- act_3.2
    di <- rbind(di,d1)
  }
  
return(di)}

# eee <- 320
# tms_3.2 <- all_tms_ex(eee)
# 
# 
# x <- c("Actual_all","CF_FBA_all","CF_FVA_min_all","CF_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all","CFMTR_FBA", "CFMTR_FVA_min","CFMTR_FVA_max" )
# si <- which(colnames(tms_3.2) %in% x)
# 
# 
# for(i in 1:length(si)){
#   tms_3.2[[si[i]]] <- scale(tms_3.2[[si[i]]], center = min(tms_3.2[[si[i]]]), scale = max(tms_3.2[[si[i]]]) - min(tms_3.2[[si[i]]]))
# }
# 
# tms_3.2_ <- tms_3.2[c(4,12,20,28,36,44,59,70,81,92,100,109,118),]
# 
# 
# library(ggplot2)
# library(ggpubr)
# 
# ggplot(tms_3.2_, aes(Actual_all, CF_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(tms_3.2_, aes(Actual_all, CFMTR_FVA_max))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(tms_3.2_, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(tms_3.2_, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# 
# 
# 
# eee <- 320
# tms_3.2 <- all_tms_ex(eee)
# 
# 
# x <- c("Actual_all","CF_FBA_all","CF_FVA_min_all","CF_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all","CFMTR_FBA", "CFMTR_FVA_min","CFMTR_FVA_max" )
# si <- which(colnames(tms_3.2) %in% x)
# 
# 
# library(caret)
# 
# tms_3.2_1 <- tms_3.2[,si]
# x <- preProcess(tms_3.2_1, method=c("range"))




############################################# ALL TMS, EXCHANGE RATES TAHEN TOGETHER






library(ggplot2)
library(ggpubr)





tms_3.2 <- all_tms_ex(3.2)
tms_320 <- all_tms_ex(320)
tms_3200 <- all_tms_ex(3200)

TMS_EXC_ALL <- rbind(tms_3.2,tms_320,tms_3200)


x <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
si <- which(colnames(TMS_EXC_ALL) %in% x)

# #
# for(i in 1:length(si)){
#   TMS_EXC_ALL[[si[i]]] <- scale(TMS_EXC_ALL[[si[i]]], center = min(TMS_EXC_ALL[[si[i]]]), scale = max(TMS_EXC_ALL[[si[i]]]) - min(TMS_EXC_ALL[[si[i]]]))
# }
# 
# 
# 
# ggplot(TMS_EXC_ALL, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(TMS_EXC_ALL, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(TMS_EXC_ALL, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
# ggplot(TMS_EXC_ALL, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various TMS under ",eee," exchange rate"))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)


############## Extracting from TMSE_EXC_ALL
ji <- c(4,12,20,28,36,44,59,70,81,92,100,109,118,
        124,131,139,147,155,163,178,189,200,211,219,228,237,
        242,250,258,266,274,282,297,308,319,330,338,347,356)
TMS_EXC_ALL_ <- TMS_EXC_ALL[ji, ]


x <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
si <- which(colnames(TMS_EXC_ALL_) %in% x)


for(i in 1:length(si)){
  TMS_EXC_ALL_[[si[i]]] <- scale(TMS_EXC_ALL_[[si[i]]], center = min(TMS_EXC_ALL_[[si[i]]]), scale = max(TMS_EXC_ALL_[[si[i]]]) - min(TMS_EXC_ALL_[[si[i]]]))
}



somed <- ("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/")


dir.create(paste0(somed,"/BM_ACT_VS_PRED/"))


setwd(paste0(somed,"/BM_ACT_VS_PRED/"))

ggplot(TMS_EXC_ALL_, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("CF_S_FVA_MAX.pdf")
ggsave("CF_S_FVA_MAX.jpeg")

ggplot(TMS_EXC_ALL_, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("CF_MTR_FVA_MAX.pdf")
ggsave("CF_MTR_FVA_MAX.jpeg")

ggplot(TMS_EXC_ALL_, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("TRIMER.pdf")
ggsave("TRIMER.jpeg")

ggplot(TMS_EXC_ALL_, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("GIMME_FVA_max.pdf")
ggsave("GIMME_FVA_max.jpeg")

# +geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
#######################################################################################################
#######################################################################################################
####### Alternate plotting  2
#######################################################################################################


t <- 3

e <- c(3.2, 320, 3200)


di <- data.frame()
for(i in 1:length(e)){
  setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TM",t,"/"))
  act_3.2 <- read.csv(paste0("Actual_Pred_data_",e[i],"_TM_",t,".csv"), header = TRUE)
  
  tm <- paste0("TM",t)
  p <- "Per_0"
  ex <- e[i]
  
  if(t == 2){
    kv <- c("WT","I","X","A")
  }else if(t == 1){
    kv <- c("WT","B","E","Z","A","X")
  }else if(t == 3){
    kv <- c("WT","A","X")
  }
  
  
  
  cfmtr_tm1_3.2_per0 <- get_overall_df(tm,p,ex,kv)
  
  act_3.2$CFMTR_FBA <- cfmtr_tm1_3.2_per0[[1]]
  act_3.2$CFMTR_FVA_min <- cfmtr_tm1_3.2_per0[[2]]
  act_3.2$CFMTR_FVA_max <- cfmtr_tm1_3.2_per0[[3]]
  
  d1 <- act_3.2
  di <- rbind(di,d1)
}

x <- c("Actual_all","CF_FBA_all","CF_FVA_min_all","CF_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all","CFMTR_FBA", "CFMTR_FVA_min","CFMTR_FVA_max" )
si <- which(colnames(di) %in% x)


for(i in 1:length(si)){
  di[[si[i]]] <- scale(di[[si[i]]], center = min(di[[si[i]]]), scale = max(di[[si[i]]]) - min(di[[si[i]]]))
}


if(t == 1){
  dii <- di[c(4,12,20,28,36,44,52,60,68,76,84,92,100,108,116,124,132,140),]
  #dii <- di[c(4,12,20,28,36,44),]
}else if(t == 2){
  dii <- di[c(11,22,33,44,55,66,77,88,99,110,121,132),]
  #dii <- di[c(11,22,33,44),]
}else if(t == 3){
  dii <- di[c(8,17,26,35,44,53,62,71,80),]
}


ggplot(dii, aes(Actual_all, CF_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, CFMTR_FVA_max))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)



#######################################################################################################
#######################################################################################################
####### Alternate plotting  3
#######################################################################################################


t <- 2

e <- c(3.2, 320, 3200)


di <- data.frame()
for(i in 1:length(e)){
  setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TM",t,"/"))
  act_3.2 <- read.csv(paste0("Actual_Pred_data_",e[i],"_TM_",t,".csv"), header = TRUE)
  
  tm <- paste0("TM",t)
  p <- "Per_0"
  ex <- e[i]
  
  if(t == 2){
    kv <- c("WT","I","X","A")
  }else if(t == 1){
    kv <- c("WT","B","E","Z","A","X")
  }else if(t == 3){
    kv <- c("WT","A","X")
  }
  
  
  
  cfmtr_tm1_3.2_per0 <- get_overall_df(tm,p,ex,kv)
  
  act_3.2$CFMTR_FBA <- cfmtr_tm1_3.2_per0[[1]]
  act_3.2$CFMTR_FVA_min <- cfmtr_tm1_3.2_per0[[2]]
  act_3.2$CFMTR_FVA_max <- cfmtr_tm1_3.2_per0[[3]]
  
  d1 <- act_3.2
  di <- rbind(di,d1)
}


if(t == 1){
  dii <- di[c(4,12,20,28,36,44,52,60,68,76,84,92,100,108,116,124,132,140),]
  #dii <- di[c(4,12,20,28,36,44),]
}else if(t == 2){
  dii <- di[c(11,22,33,44,55,66,77,88,99,110,121,132),]
  #dii <- di[c(11,22,33,44),]
}else if(t == 3){
  dii <- di[c(8,17,26,35,44,53,62,71,80),]
}


x <- c("Actual_all","CF_FBA_all","CF_FVA_min_all","CF_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all","CFMTR_FBA", "CFMTR_FVA_min","CFMTR_FVA_max" )
si <- which(colnames(dii) %in% x)


for(i in 1:length(si)){
  dii[[si[i]]] <- scale(dii[[si[i]]], center = min(dii[[si[i]]]), scale = max(dii[[si[i]]]) - min(dii[[si[i]]]))
}





ggplot(dii, aes(Actual_all, CF_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)
ggplot(dii, aes(Actual_all, CFMTR_FVA_max))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases and various exchange rates in TM",t))+theme_bw()+geom_jitter(width = 0.05, height = 0.05, alpha= 0.25)




###########################################################################################################################
###########################################################################################################################

#### RMSE from TMS_EXC_ALL

library(Metrics)

############# Overall 
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$CF_FVA_max_all) # 1072.141
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$CFMTR_FVA_max) # 1068.465
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$GIMME_FVA_max_all) # 1133.981
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$TRIMER_all) # 676.9728


############ with feedback 
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$CF_FVA_max_all) # 1140.553
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$CFMTR_FVA_max) # 1095.507
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$GIMME_FVA_max_all) # 1316.469
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$TRIMER_all) # 802.6443


########## overall after normalization 
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$CF_FVA_max_all) # 0.2599125
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$CFMTR_FVA_max) # 0.2638284
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$GIMME_FVA_max_all) # 0.2673023
rmse(TMS_EXC_ALL$Actual_all, TMS_EXC_ALL$TRIMER_all) # 0.1945571


############ with feedback after normalization 
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$CF_FVA_max_all) # 0.2854848
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$CFMTR_FVA_max) # 0.2739648
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$GIMME_FVA_max_all) # 0.3388376
rmse(TMS_EXC_ALL_$Actual_all, TMS_EXC_ALL_$TRIMER_all) # 0.2079335





##########
setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/Results_1/BM_act_vs_pred/new/")
write.csv(TMS_EXC_ALL, file = "TMS_EXC_ALL.csv", quote = FALSE, col.names = TRUE, row.names = FALSE)


################################################################################################################
################################################################################################################
################################################################################################################
#### New BM act vs  pred

################################################################################################################
################################################################################################################
################################################################################################################


obtain_df <- function(tm, ex, co){
  library(readxl)
  setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/",tm,"/Results/",ex,"/KO_data_",tm,"_",ex,"/",co))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_WT <- y$V1
  
  z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/",tm,"/Results/",ex,"/KO_data_",tm,"_",ex,"/",co,"/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_WT <- z[,2]
  FVA_max_WT <- z[,3]
  
  tg_df <- data.frame(FBA_WT, FVA_min_WT, FVA_max_WT)
  colnames(tg_df) <- c("FBA","FVA_min","FVA_max") 
  tg_df$condition <- rep(co, nrow(tg_df))
  
  return(tg_df)}


get_overall_df <- function(tm, ex, kv){
  ddf <- data.frame()
  
  for(i in 1:length(kv)){
    wt <- obtain_df(tm,ex,kv[i])
    ddf <- rbind(ddf,wt)
  }
  
  
  return(ddf)}



all_tms_ex <- function(e){
  t <- c(1, 2, 3)
  di <- data.frame()
  
  for(i in 1:length(t)){
    setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/TM",t[i],"/"))
    act_3.2 <- read.csv(paste0("Actual_Pred_data_",e,"_TM_",t[i],".csv"), header = TRUE)
    
    tm <- paste0("TM",t[i])
    #p <- "Per_0"
    ex <- e
    
    if(t[i] == 2){
      kv <- c("WT","I","X","A")
    }else if(t[i] == 1){
      kv <- c("WT","B","E","Z","A","X")
    }else if(t[i] == 3){
      kv <- c("WT","A","X")
    }
    
    
    
    cfmtr_tm1_3.2_per0 <- get_overall_df(tm,ex,kv)
    
    act_3.2$CFMTR_FBA <- cfmtr_tm1_3.2_per0[[1]]
    act_3.2$CFMTR_FVA_min <- cfmtr_tm1_3.2_per0[[2]]
    act_3.2$CFMTR_FVA_max <- cfmtr_tm1_3.2_per0[[3]]
    
    d1 <- act_3.2
    di <- rbind(di,d1)
  }
  
  return(di)}



tms_3.2 <- all_tms_ex(3.2)
tms_320 <- all_tms_ex(320)
tms_3200 <- all_tms_ex(3200)

TMS_EXC_ALL <- rbind(tms_3.2,tms_320,tms_3200)

ji <- c(4,12,20,28,36,44,59,70,81,92,100,109,118,
        124,131,139,147,155,163,178,189,200,211,219,228,237,
        242,250,258,266,274,282,297,308,319,330,338,347,356)
TMS_EXC_ALL_ <- TMS_EXC_ALL[ji, ]



ggplot(TMS_EXC_ALL_, aes(Actual_all, CF_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()
ggplot(TMS_EXC_ALL_, aes(Actual_all, CFMTR_FVA_max))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()
ggplot(TMS_EXC_ALL_, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()
ggplot(TMS_EXC_ALL_, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()









