###############################################################################################
###############################################################################################


curr_wd <- c("D:/work/Integrated_network_model/Git_hub_codes/Ecoli/PROM_rFBA_CausalFlux/")

#############################################################################################
####################  PROM  ---- 0% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] == 0*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"PROM/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] ==0 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] ==0 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)



ggplot(tot_comp_df, aes(x = cfs_act, y = act_prom)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs PROM (Overall media) - 0% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between PROM and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_prom_0per_wt.pdf")
ggsave("cf_prom_0per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

# tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
# otp <- c(otp,tot_cm_cfs$table[1,1])
# ofp <- c(ofp,tot_cm_cfs$table[1,2])
# ofn <- c(ofn,tot_cm_cfs$table[2,1])
# otn <- c(otn,tot_cm_cfs$table[2,2])
# op <- c(op,tot_cm_cfs$byClass["Precision"])
# or <- c(or, tot_cm_cfs$byClass["Recall"])
# of <- c(of,  tot_cm_cfs$byClass["F1"])
# oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])


tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_prom_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf1 <- data.frame(
  "Methods" = c("CausalFlux pred on PROM dataset (0% of WT)","PROM"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)


#############################################################################################
####################  PROM  ---- 5% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"PROM/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)


ggplot(tot_comp_df, aes(x = cfs_act, y = act_prom)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs PROM (Overall media) - 5% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between PROM and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_prom_5per_wt.pdf")
ggsave("cf_prom_5per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

# tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
# otp <- c(otp,tot_cm_cfs$table[1,1])
# ofp <- c(ofp,tot_cm_cfs$table[1,2])
# ofn <- c(ofn,tot_cm_cfs$table[2,1])
# otn <- c(otn,tot_cm_cfs$table[2,2])
# op <- c(op,tot_cm_cfs$byClass["Precision"])
# or <- c(or, tot_cm_cfs$byClass["Recall"])
# of <- c(of,  tot_cm_cfs$byClass["F1"])
# oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])


tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_prom_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf2 <- data.frame(
  "Methods" = c("CausalFlux pred on PROM dataset (5% of WT)","PROM"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)


#############################################################################################
####################  PROM  ---- 50% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"PROM/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"PROM/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"PROM/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)


ggplot(tot_comp_df, aes(x = cfs_act, y = act_prom)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs PROM (Overall media) - 50% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between PROM and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_prom_50per_wt.pdf")
ggsave("cf_prom_50per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

# tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
# otp <- c(otp,tot_cm_cfs$table[1,1])
# ofp <- c(ofp,tot_cm_cfs$table[1,2])
# ofn <- c(ofn,tot_cm_cfs$table[2,1])
# otn <- c(otn,tot_cm_cfs$table[2,2])
# op <- c(op,tot_cm_cfs$byClass["Precision"])
# or <- c(or, tot_cm_cfs$byClass["Recall"])
# of <- c(of,  tot_cm_cfs$byClass["F1"])
# oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])


tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_prom_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf3 <- data.frame(
  "Methods" = c("CausalFlux pred on PROM dataset (50% of WT)","PROM"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)



#############################################################################################
### rFBA
#############################################################################################

#############################################################################################
####################  rFBA  ---- 0% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] == 0*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"rFBA/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] ==0 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] ==0 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)


setwd(curr_wd)
ggplot(tot_comp_df, aes(x = cfs_act, y = act_rfba)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs rFBA (Overall media) - 0% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between rFBA and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_rfba_0per_wt.pdf")
ggsave("cf_rfba_0per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf4 <- data.frame(
  "Methods" = c("CausalFlux pred on rFBA dataset (0% of WT)","rFBA"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)


#############################################################################################
####################  rFBA  ---- 5% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"rFBA/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.05 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)


setwd(curr_wd)
ggplot(tot_comp_df, aes(x = cfs_act, y = act_rfba)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs rFBA (Overall media) - 5% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between rFBA and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_rfba_5per_wt.pdf")
ggsave("cf_rfba_5per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf5 <- data.frame(
  "Methods" = c("CausalFlux pred on rFBA dataset (5% of WT)","rFBA"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)

#############################################################################################
####################  rFBA  ---- 50% of WT 

######################### Carbon Source
med_con_1 <- c("CF_S_rfba_EX_12ppd__S_e","CF_S_rfba_EX_dad_2_e","CF_S_rfba_EX_glc__D_e","CF_S_rfba_EX_lcts_e","CF_S_rfba_EX_akg_e",
               "CF_S_rfba_EX_ac_e",
               "CF_S_rfba_EX_acac_e","CF_S_rfba_EX_adn_e","CF_S_rfba_EX_cit_e","CF_S_rfba_EX_mal__D_e","CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_fru_e","CF_S_rfba_EX_gal_e","CF_S_rfba_EX_galur_e","CF_S_rfba_EX_glcn_e",
               "CF_S_rfba_EX_g6p_e",
               "CF_S_rfba_EX_mnl_e","CF_S_rfba_EX_man_e","CF_S_rfba_EX_melib_e",
               "CF_S_rfba_EX_rib__D_e",
               "CF_S_rfba_EX_ser__D_e","CF_S_rfba_EX_sbt__D_e","CF_S_rfba_EX_tre_e",
               "CF_S_rfba_EX_xyl__D_e",
               "CF_S_rfba_EX_for_e",
               "CF_S_rfba_EX_fum_e","CF_S_rfba_EX_glyc_e",
               "CF_S_rfba_EX_glyclt_e",
               "CF_S_rfba_EX_ins_e","CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_arab__L_e","CF_S_rfba_EX_asn__L_e","CF_S_rfba_EX_asp__L_e",
               "CF_S_rfba_EX_glu__L_e",
               "CF_S_rfba_EX_gln__L_e",
               "CF_S_rfba_EX_lac__L_e","CF_S_rfba_EX_mal__L_e","CF_S_rfba_EX_pro__L_e",
               "CF_S_rfba_EX_rmn_e","CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e","CF_S_rfba_EX_malt_e","CF_S_rfba_EX_malttr_e","CF_S_rfba_EX_acmana_e","CF_S_rfba_EX_acgam_e", 
               "CF_S_rfba_EX_pyr_e",
               "CF_S_rfba_EX_succ_e","CF_S_rfba_EX_sucr_e","CF_S_rfba_EX_thymd_e","CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_4abut_e", "CF_S_rfba_EX_acnam_e",
               "CF_S_rfba_EX_arg__L_e","CF_S_rfba_EX_but_e","CF_S_rfba_EX_dha_e",
               "CF_S_rfba_EX_orn_e","CF_S_rfba_EX_ptrc_e","CF_S_rfba_EX_tartr__L_e")



gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_CS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))

for(j in 1:length(med_con_1)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    #new_Iterations_10
    setwd(paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/CS/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    
    xx_end <- max(file_numbers)
    
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    # 
    
    # val_r_1 <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/Separate_PROM_0p/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_1b_obj_0_P1.xlsx"), col_names = F)  
    # xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_CS <- rbind(run_set_CS,med_1)
}


PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5*wtb){
        #
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805

PN_df_chk <- PN_df(run_set_CS,WT_bm)

setwd(paste0(curr_wd,"rFBA/CS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_cs <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_cs <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_cs <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk, df_actual_cs)

cfs_rfba <- diff_per_cal(PN_df_chk, df_rfba_cs)
cfs_prom <- diff_per_cal(PN_df_chk, df_prom_cs)

act_prom <- diff_per_cal(df_actual_cs, df_prom_cs)
act_rfba <- diff_per_cal(df_actual_cs, df_rfba_cs)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set$medconds <- overall_GR$Growth.Media


colnames(run_set_CS) <- gene_ko
pred_cfs <- as.vector(t(run_set_CS))

cv_cfs <- as.vector(t(PN_df_chk))
cv_actual <- as.vector(t(df_actual_cs))
cv_rfba <- as.vector(t(df_rfba_cs))
cv_prom <- as.vector(t(df_prom_cs))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Carbon_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)

#########################  Nitrogen source


med_con_2 <- c("CF_S_rfba_EX_ade_e",
               "CF_S_rfba_EX_adn_e",
               "CF_S_rfba_EX_alltn_e",
               "CF_S_rfba_EX_cytd_e",
               "CF_S_rfba_EX_ala__D_e",
               "CF_S_rfba_EX_gam_e",
               "CF_S_rfba_EX_gsn_e",
               "CF_S_rfba_EX_ins_e",
               "CF_S_rfba_EX_ala__L_e",
               "CF_S_rfba_EX_leu__L_e",
               "CF_S_rfba_EX_lys__L_e",
               "CF_S_rfba_EX_met__D_e",
               "CF_S_rfba_EX_ser__L_e",
               "CF_S_rfba_EX_thr__L_e",
               "CF_S_rfba_EX_trp__L_e",
               "CF_S_rfba_EX_tyr__L_e",
               "CF_S_rfba_EX_ptrc_e",
               "CF_S_rfba_EX_thymd_e",
               "CF_S_rfba_EX_ura_e",
               "CF_S_rfba_EX_urea_e",
               "CF_S_rfba_EX_uri_e",
               "CF_S_rfba_EX_xan_e",
               "CF_S_rfba_EX_xtsn_e",
               "CF_S_rfba_EX_acmana_e",
               "CF_S_rfba_EX_acgam_e",
               "CF_S_rfba_EX_gua_e",
               "CF_S_rfba_EX_ile__L_e",
               "CF_S_rfba_EX_phe__L_e",
               "CF_S_rfba_EX_gly_e",
               "CF_S_rfba_EX_cys__L_e",
               "CF_S_rfba_EX_his__L_e"
               
)


gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko


run_set_NS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))


for(j in 1:length(med_con_2)){ 
  med_2 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/NS/",med_con_2[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_2 <- c(med_2,((xx_val_end+xx_val_end_)/2))
    
    med_2 <- c(med_2,xx_val_end)
  }
  run_set_NS <- rbind(run_set_NS,med_2)
}

PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5 * wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_ns <- PN_df(run_set_NS,WT_bm)

############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/NS/"))
overall_GR <- read.csv("run_set_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_ns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_ns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_ns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_ns, df_actual_ns)

cfs_rfba <- diff_per_cal(PN_df_chk_ns, df_rfba_ns)
cfs_prom <- diff_per_cal(PN_df_chk_ns, df_prom_ns)

act_prom <- diff_per_cal(df_actual_ns, df_prom_ns)
act_rfba <- diff_per_cal(df_actual_ns, df_rfba_ns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_ns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_ns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_NS))

cv_cfs <- as.vector(t(PN_df_chk_ns))
cv_actual <- as.vector(t(df_actual_ns))
cv_rfba <- as.vector(t(df_rfba_ns))
cv_prom <- as.vector(t(df_prom_ns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)


################################### Double Nitrogen source

med_con_3 <- c("CF_S_rfba_EX_ala_EX_asp__L_e","CF_S_rfba_EX_ala_EX_gln__L_e","CF_S_rfba_EX_ala_EX_glu__L_e",
               "CF_S_rfba_EX_ala_EX_gly_e", "CF_S_rfba_EX_ala_EX_his__L_e", "CF_S_rfba_EX_ala_EX_leu__L_e",
               "CF_S_rfba_EX_ala_EX_thr__L_e", "CF_S_rfba_EX_gly_EX_gln__L_e", "CF_S_rfba_EX_gly_EX_glu__L_e","CF_S_rfba_EX_gly_EX_met__L_e" ,"CF_S_rfba_EX_met_EX_ala__L_e")




gene_ko <- c("tdcR","crp", "malT","glpR","gntR","xylR","asnC","rbsR","ilvY","glnG","rhaS","cpxR","cytR","soxR","melR")

m <- gene_ko

##################################### med cond 1
run_set_DNS <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(gene_ko)))



for(j in 1:length(med_con_3)){ 
  med_1 <- c()
  for(i in 1:length(gene_ko)){
    #val_r <- readxl::read_xlsx(paste0("D:/work/Integrated_network_model/papers_materials/rFBA/Runs/Run_Results/new_results_based_on_prom_grn/modified_upper_bounds_Carbon_Source/run_Set_2/",med_con_1[j],"/Case_P1_",i,"/",m[i],"/FVA_to_check_P1.xlsx"), col_names = F)
    
    setwd(paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/"))
    
    folder_path <- paste0(curr_wd,"rFBA/DNS/",med_con_3[j],"/Case_P1_",i,"/",m[i],"/")
    
    # List all csv files starting with 'F_'
    files <- list.files(folder_path, pattern = "^FVA_incorp_P1_\\d+\\.csv$", full.names = TRUE)
    
    # Extract the numeric part (e.g., 1, 2, 3, ...) using regex
    file_numbers <- as.numeric(gsub(".*FVA_incorp_P1_(\\d+)\\.csv", "\\1", files))
    
    
    xx_end <- max(file_numbers)
    #xx_end <- 1
    
    val_r_1 <- read.csv(paste0("FVA_incorp_P1_",xx_end,".csv"))
    xx_val_end <- val_r_1[[3]][269]
    
    # xx_end_ <- max(file_numbers)-1
    # val_r_2 <- read.csv(paste0("FVA_incorp_P1_",xx_end_,".csv"))
    # xx_val_end_ <- val_r_2[[3]][269]
    
    
    #med_1 <- c(med_1,((xx_val_end+xx_val_end_)/2))
    
    med_1 <- c(med_1,xx_val_end)
  }
  run_set_DNS <- rbind(run_set_DNS,med_1)
}



PN_df <- function(ddf,wtb){
  for(i in 1:nrow(ddf)){
    for(j in 1:ncol(ddf)){
      if(ddf[i,j] <= 0.5 *wtb){
        ddf[i,j] <- "N"
      }else{
        ddf[i,j] <- "P"
      }
    }
  }
  return(ddf)}

WT_bm <- 0.585503805
PN_df_chk_dns <- PN_df(run_set_DNS,WT_bm)


############################# read the actual, prom, rfba data 
setwd(paste0(curr_wd,"rFBA/DNS/"))
overall_GR <- read.csv("NS_2.csv", header = T)

OGR <- overall_GR[,-1]

OGR_new <- as.data.frame(apply(OGR, 2, function(x) {
  ifelse(x == "+", "P", ifelse(x == "-", "N", x))
}), stringsAsFactors = FALSE)

# Columns ending with "_a"
df_actual_dns <- OGR_new[ , grep("_actual$", colnames(OGR_new)) ]

# Columns ending with "_b"
df_prom_dns <- OGR_new[ , grep("_prom$", colnames(OGR_new)) ]

# Columns ending with "_c"
df_rfba_dns <- OGR_new[ , grep("_rfba$", colnames(OGR_new)) ]



diff_per_cal <- function(d1, d2){
  dfp <-c()
  for(i in 1:nrow(d1)){
    dfp <- c(dfp,mean(d1[i,]==d2[i,]) * 100)
  }
  
  return(dfp)}

cfs_act <- diff_per_cal(PN_df_chk_dns, df_actual_dns)

cfs_rfba <- diff_per_cal(PN_df_chk_dns, df_rfba_dns)
cfs_prom <- diff_per_cal(PN_df_chk_dns, df_prom_dns)

act_prom <- diff_per_cal(df_actual_dns, df_prom_dns)
act_rfba <- diff_per_cal(df_actual_dns, df_rfba_dns)
#prom_rfba <- diff_per_cal(prom_pred, rfba_pred)

comp_df_run_set_dns <- data.frame(cbind(cfs_act,cfs_rfba,cfs_prom,act_prom,act_rfba))

comp_df_run_set_dns$medconds <- overall_GR$Growth.Media



pred_cfs <- as.vector(t(run_set_DNS))

cv_cfs <- as.vector(t(PN_df_chk_dns))
cv_actual <- as.vector(t(df_actual_dns))
cv_rfba <- as.vector(t(df_rfba_dns))
cv_prom <- as.vector(t(df_prom_dns))


cv_cfs_nv <- ifelse(cv_cfs == "P", 1, 0)
cv_actual_nv <- ifelse(cv_actual == "P", 1, 0)
cv_rfba_nv <- ifelse(cv_rfba == "P", 1, 0)
cv_prom_nv <- ifelse(cv_prom == "P", 1, 0)



Double_Nitrogen_Souce_DF <- data.frame(cv_cfs_nv,cv_actual_nv,cv_rfba_nv,cv_prom_nv)





tot_comp_df <- rbind(comp_df_run_set, comp_df_run_set_ns, comp_df_run_set_dns)


ggplot(tot_comp_df, aes(x = cfs_act, y = act_rfba)) +
  geom_abline() +
  geom_jitter(width = 1.75, height = 1.75, alpha = 0.4, color = "firebrick", size = 3) +
  ggtitle("Jitter plot: CausalFlux vs rFBA (Overall media) - 50% of WT") +
  xlab("Accuracy between CausalFlux and Actual") +
  ylab("Accuracy between rFBA and Actual") +
  theme_bw()

setwd(curr_wd)
ggsave("cf_rfba_50per_wt.pdf")
ggsave("cf_rfba_50per_wt.jpeg")



Tot_df <- rbind(Carbon_Souce_DF,Nitrogen_Souce_DF,Double_Nitrogen_Souce_DF)
otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_cfs_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])

tot_cm_cfs <- confusionMatrix(as.factor(Tot_df$cv_rfba_nv), as.factor(Tot_df$cv_actual_nv), positive = "0")
otp <- c(otp,tot_cm_cfs$table[1,1])
ofp <- c(ofp,tot_cm_cfs$table[1,2])
ofn <- c(ofn,tot_cm_cfs$table[2,1])
otn <- c(otn,tot_cm_cfs$table[2,2])
op <- c(op,tot_cm_cfs$byClass["Precision"])
or <- c(or, tot_cm_cfs$byClass["Recall"])
of <- c(of,  tot_cm_cfs$byClass["F1"])
oba <- c(oba,  tot_cm_cfs$byClass["Balanced Accuracy"])



somedf6 <- data.frame(
  "Methods" = c("CausalFlux pred on rFBA dataset (50% of WT)","rFBA"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)





odf <- rbind(somedf1[1,], somedf2[1,], somedf3[1,],somedf4[1,],somedf5[1,],somedf6[1,], somedf1[2,], somedf4[2,])
setwd(curr_wd)
write.csv(odf,"Metrics_CausalFlux_PROM_rFBA.csv")



