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



