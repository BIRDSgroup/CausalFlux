######################################################################################################
######################  Functions for this code
get_overall_df <- function(tm, ex, kv){
  ddf <- data.frame()
  
  for(i in 1:length(kv)){
    wt <- obtain_df(tm,ex,kv[i])
    ddf <- rbind(ddf,wt)
  }
  
  
  return(ddf)}


#############################################################################################
##############  TM1 

act_vs_pred_sep_plots_tm1 <- function(ec,g, somed){
  
  e <- c(3.2,320,3200)
  
  if(ec == 3.2){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_1.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_1.csv", header = TRUE)

    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_1.csv", header = TRUE)
    
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
  #x2 <- c("Actual","CF-S (FBA)","CF-S (FVA min)","CF-S (FVA max)", "TRIMER","GIMME (FBA)","GIMME (FVA min)","GIMME (FVA max)","CF-MTR (FBA)","CF-MTR (FVA min)","CF-MTR (FVA max)")
  
  #x <- data.frame(x1,x2)
  #si <- which(colnames(d1) %in% x[,1])
  si <- which(colnames(d1) %in% x1)

  for(i in 1:length(si)){
    d1[[si[i]]] <- scale(d1[[si[i]]], center = min(d1[[si[i]]]), scale = max(d1[[si[i]]]) - min(d1[[si[i]]]))
  }

  
  
  library(ff)
  
  #somed <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Results/Sep_act_vs_pred_plot/"
  
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot","/TM1"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM1/", ec))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM1/", ec, "/",g))
  
  EO <- ec
  EM <- ec
  COND <- g
  TM = "TM 1"
  
  library(ggplot2)
  library(ggpubr)
  
  setwd(paste0(somed,"/Sep_act_vs_pred_plot/TM1/", ec, "/",g))
  
  
  ggplot(d1, aes(Actual_all, CF_S_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FBA.pdf")
  ggsave("CF-S_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MIN.pdf")
  ggsave("CF-S_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MAX.pdf")
  ggsave("CF-S_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("TRIMER.pdf")
  ggsave("TRIMER.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FBA.pdf")
  ggsave("GIMME_FBA.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MIN.pdf")
  ggsave("GIMME_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MAX.pdf")
  ggsave("GIMME_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FBA.pdf")
  ggsave("CF_MTR_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MIN.pdf")
  ggsave("CF_MTR_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MAX.pdf")
  ggsave("CF_MTR_FVA_MAX.jpeg")
  
  
  
}


############################################# 3.2
#sd <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Results/Sep_act_vs_pred_plot/"

sd <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/"



act_vs_pred_sep_plots_tm1(3.2,"WT",sd)
act_vs_pred_sep_plots_tm1(3.2,"X KO",sd)
act_vs_pred_sep_plots_tm1(3.2,"A KO",sd)
act_vs_pred_sep_plots_tm1(3.2,"B KO",sd)
act_vs_pred_sep_plots_tm1(3.2,"E KO",sd)
act_vs_pred_sep_plots_tm1(3.2,"Z KO",sd)


############################################# 320
act_vs_pred_sep_plots_tm1(320,"WT",sd)
act_vs_pred_sep_plots_tm1(320,"X KO",sd)
act_vs_pred_sep_plots_tm1(320,"A KO",sd)
act_vs_pred_sep_plots_tm1(320,"B KO",sd)
act_vs_pred_sep_plots_tm1(320,"E KO",sd)
act_vs_pred_sep_plots_tm1(320,"Z KO",sd)

############################################# 3200
act_vs_pred_sep_plots_tm1(3200,"WT",sd)
act_vs_pred_sep_plots_tm1(3200,"X KO",sd)
act_vs_pred_sep_plots_tm1(3200,"A KO",sd)
act_vs_pred_sep_plots_tm1(3200,"B KO",sd)
act_vs_pred_sep_plots_tm1(3200,"E KO",sd)
act_vs_pred_sep_plots_tm1(3200,"Z KO",sd)



#############################################################################################
##############  TM2 

act_vs_pred_sep_plots_tm2 <- function(ec,g, somed){
  
  e <- c(3.2,320,3200)
  
  if(ec == 3.2){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
  #x2 <- c("Actual","CF-S (FBA)","CF-S (FVA min)","CF-S (FVA max)", "TRIMER","GIMME (FBA)","GIMME (FVA min)","GIMME (FVA max)","CF-MTR (FBA)","CF-MTR (FVA min)","CF-MTR (FVA max)")
  
  #x <- data.frame(x1,x2)
  #si <- which(colnames(d1) %in% x[,1])
  si <- which(colnames(d1) %in% x1)
  
  for(i in 1:length(si)){
    d1[[si[i]]] <- scale(d1[[si[i]]], center = min(d1[[si[i]]]), scale = max(d1[[si[i]]]) - min(d1[[si[i]]]))
  }
  
  
  
  library(ff)
  
  #somed <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/Results_1/Sep_act_vs_pred_plot"
  
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot","/TM2"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM2/", ec))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM2/", ec, "/",g))
  
  EO <- ec
  EM <- ec
  COND <- g
  TM = "TM 2"
  
  library(ggplot2)
  library(ggpubr)
  
  setwd(paste0(somed,"/Sep_act_vs_pred_plot/TM2/", ec, "/",g))
  
  
  ggplot(d1, aes(Actual_all, CF_S_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FBA.pdf")
  ggsave("CF-S_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MIN.pdf")
  ggsave("CF-S_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MAX.pdf")
  ggsave("CF-S_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("TRIMER.pdf")
  ggsave("TRIMER.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FBA.pdf")
  ggsave("GIMME_FBA.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MIN.pdf")
  ggsave("GIMME_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MAX.pdf")
  ggsave("GIMME_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FBA.pdf")
  ggsave("CF_MTR_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MIN.pdf")
  ggsave("CF_MTR_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MAX.pdf")
  ggsave("CF_MTR_FVA_MAX.jpeg")
  
  
  
}

############################################# 320
act_vs_pred_sep_plots_tm2(320,"X KO",sd)

act_vs_pred_sep_plots_tm2(320,"WT",sd)

act_vs_pred_sep_plots_tm2(320,"A KO",sd)

act_vs_pred_sep_plots_tm2(320,"I KO",sd)


############################################# 3.2
act_vs_pred_sep_plots_tm2(3.2,"X KO",sd)

act_vs_pred_sep_plots_tm2(3.2,"WT",sd)

act_vs_pred_sep_plots_tm2(3.2,"A KO",sd)

act_vs_pred_sep_plots_tm2(3.2,"I KO",sd)


############################################# 3200
act_vs_pred_sep_plots_tm2(3200,"X KO",sd)

act_vs_pred_sep_plots_tm2(3200,"WT",sd)

act_vs_pred_sep_plots_tm2(3200,"A KO",sd)

act_vs_pred_sep_plots_tm2(3200,"I KO",sd)




#############################################################################################
##############  TM3 

act_vs_pred_sep_plots_tm3 <- function(ec,g, somed){
  
  e <- c(3.2,320,3200)
  
  if(ec == 3.2){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Final_simulation_results/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","CF_M_FBA_all", "CF_M_FVA_min_all","CF_M_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
  #x2 <- c("Actual","CF-S (FBA)","CF-S (FVA min)","CF-S (FVA max)", "TRIMER","GIMME (FBA)","GIMME (FVA min)","GIMME (FVA max)","CF-MTR (FBA)","CF-MTR (FVA min)","CF-MTR (FVA max)")
  
  #x <- data.frame(x1,x2)
  #si <- which(colnames(d1) %in% x[,1])
  si <- which(colnames(d1) %in% x1)
  
  for(i in 1:length(si)){
    d1[[si[i]]] <- scale(d1[[si[i]]], center = min(d1[[si[i]]]), scale = max(d1[[si[i]]]) - min(d1[[si[i]]]))
  }
  
  
  
  library(ff)
  
  #somed <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/CF_MTR/Results_1/Sep_act_vs_pred_plot/"
  
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot","/TM3"))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM3/", ec))
  dir.create(paste0(somed,"/Sep_act_vs_pred_plot/TM3/", ec, "/",g))
  
  EO <- ec
  EM <- ec
  COND <- g
  TM = "TM 3"
  
  library(ggplot2)
  library(ggpubr)
  
  setwd(paste0(somed,"/Sep_act_vs_pred_plot/TM3/", ec, "/",g))
  
  
  ggplot(d1, aes(Actual_all, CF_S_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FBA.pdf")
  ggsave("CF-S_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MIN.pdf")
  ggsave("CF-S_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF-S_FVA_MAX.pdf")
  ggsave("CF-S_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("TRIMER.pdf")
  ggsave("TRIMER.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FBA.pdf")
  ggsave("GIMME_FBA.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MIN.pdf")
  ggsave("GIMME_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("GIMME_FVA_MAX.pdf")
  ggsave("GIMME_FVA_MAX.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FBA.pdf")
  ggsave("CF_MTR_FBA.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MIN.pdf")
  ggsave("CF_MTR_FVA_MIN.jpeg")
  ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  ggsave("CF_MTR_FVA_MAX.pdf")
  ggsave("CF_MTR_FVA_MAX.jpeg")
  
  
  
}


############################################# 320
act_vs_pred_sep_plots_tm3(320,"X KO",sd)

act_vs_pred_sep_plots_tm3(320,"WT",sd)

act_vs_pred_sep_plots_tm3(320,"A KO",sd)


############################################# 3.2
act_vs_pred_sep_plots_tm3(3.2,"X KO",sd)

act_vs_pred_sep_plots_tm3(3.2,"WT",sd)

act_vs_pred_sep_plots_tm3(3.2,"A KO",sd)

############################################# 3200
act_vs_pred_sep_plots_tm3(3200,"X KO",sd)

act_vs_pred_sep_plots_tm3(3200,"WT",sd)

act_vs_pred_sep_plots_tm3(3200,"A KO",sd)


















