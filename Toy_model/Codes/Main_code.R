########################################################################################
########## Automatic script ############################################################


##### For the toy model (TM) ############################################################

### Include the required functions here


######################### Step 1a -  GRN structure and parameter learning
## function to create 6 random variables with different probabilities for their data distribution; find their condition probabilities (CP);
## final output is a table with genes(nodes) and their CPs

GE_bin_df  <- function(Gene_name, Prob_initial, numb){
  
  mat <- matrix(nrow = numb, ncol = length(Gene_name))
  ge_data <- data.frame(mat)
  colnames(ge_data) <- Gene_name
  x = length(Prob_initial)
  for(i in 1:x){
    ge_data[i] <- rbinom(n=numb, size = 1, prob = Prob_initial[i])
  }
  base::return(ge_data)}


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



############################################# Have all the functions set up here
###########   Step - 2b  Binarization of flux (max from FVA) vector 
# create a function to get the binary data

# inputs to this function are: met_gene_reg_arcA and FVA_arcA_r1

get_the_binary_data <- function(df1, df2){
  bin_vec <- rep(0,times=nrow(df1))
  for(i in 1:nrow(df1)){
    for(j in 1:nrow(df2)){
      if(df1[[5]][i]==df2[[1]][j] && df2[[3]][j]>0){
        bin_vec[i]=1
      }
    }
  }
  df1$bin_data <- bin_vec
  return(df1)}


####################################### Step - 2c  Integration into GRN  ####################################### 

######################################################## This section is divided into two functions 

### Function 1 - to output the functional TF regulation
Integration_part_1 <- function(df_1, bn_obj, GS){
  
  # removing the unkown interactions
  q_idx <- c()
  for(i in 1:nrow(df_1)){
    if(df_1[[4]][i]=="?"){
      q_idx <- c(q_idx,i)
    }
  }
  
  if(length(q_idx)!=0){
    df_1 <- df_1[-q_idx,]
  }else{
    df_1 <- df_1
  }
  
  
  # to remove those TFs with both "+" and "-" interactions
  # get the unique TFs
  unique_tf <- unique(df_1[[3]])
  
  # use the TFs to figure out the indices of where they are in df_1
  UNI_TF_name_ids <- list()
  for(j in 1:length(unique_tf)){
    UNI_TF_name_ids[[j]] <- which(df_1[[3]]==unique_tf[j]) 
  }
  
  # name the list with TFs
  names(UNI_TF_name_ids) <- unique_tf
  
  # To check for those TFs that have more than one index
  TF_len_more_1_indx <- c()
  x_id <- c()
  for(i in 1:length(unique_tf)){
    if(length(UNI_TF_name_ids[[i]])>1){
      x_id <- c(x_id,i)
    }
  }
  
  more_1_gene <- c()
  for(i in 1:length(x_id)){
    more_1_gene <- c(more_1_gene,names(UNI_TF_name_ids[x_id[i]]))
  }
  
  more_1_gene_int <- list()
  for(i in 1:length(x_id)){
    more_1_gene_int[[i]] <- df_1[[4]][UNI_TF_name_ids[[x_id[i]]]]
  }
  
  names(more_1_gene_int) <- more_1_gene
  #some_ls <- list(x_id,more_1_gene_int)
  #names(some_ls) <- c("len_more_1_ids","list")
  
  # to get a df with only 1 rep of regulation
  only_1_idx <- c()
  for(i in 1:length(UNI_TF_name_ids)){
    if(length(UNI_TF_name_ids[[i]])==1){
      only_1_idx <- c(only_1_idx,i)
    }
  }
  
  req_names <- names(UNI_TF_name_ids)[only_1_idx]
  
  x_int <- df_1[which(df_1[[3]] %in% req_names),3:4]
  
  
  
  
  to_filter_reg <- function(s_list){
    
    more_2_uni <- c()
    no_2 <- c()
    no_2_reg <- c()
    for(i in 1:length(s_list)){
      if(length(unique(s_list[[i]]))>1){
        more_2_uni <- c(more_2_uni,names(s_list)[i])
      } else{
        no_2 <- c(no_2,names(s_list)[i])
        no_2_reg <- c(no_2_reg,s_list[[i]][1])
      }
    }
    no_df <- data.frame(no_2, no_2_reg)
    colnames(no_df) <- c("TF-gene","Interaction")
    
    op_list <- list(more_2_uni, no_df)
    names(op_list) <- c("TFs_to_be_removed","Multiple_reg_df")
    return(op_list)}
  
  remaining_regs <- to_filter_reg(more_1_gene_int)
  x_int_ <- remaining_regs$Multiple_reg_df
  
  
  
  
  # have a feeling this is not required
  tf_0_bv <- c()
  for(i in 1:nrow(df_1)){
    if(df_1[[6]][i]==0){
      tf_0_bv <- c(tf_0_bv, df_1[[3]][i]) 
    }
  }
  
  
  p <- which(x_int[[1]] %in% tf_0_bv)
  
  if(length(p)!=0){
    x_int <- x_int[-p,]
  }else{
    x_int <- x_int
  }
  
  
  modified_interaction_df <- rbind(x_int, x_int_)
  
  bn_arcs <- bn_obj[[3]]
  bn_arcs <- as.data.frame(bn_arcs)
  
  
  
  ## Get the target gene indices
  target_genes_list_indices <- list()
  for(i in 1:length(modified_interaction_df[[1]])){
    for(j in 1:length(bn_arcs$from)){
      target_genes_list_indices[[i]] <- which(bn_arcs$from==modified_interaction_df[[1]][i])
    }
  }
  names(target_genes_list_indices) <- modified_interaction_df[[1]]
  
  ## Get the target gene names
  
  target_genes_list_names <- list()
  for(i in 1:length(target_genes_list_indices)){
    target_genes_list_names[[i]] <- as.character(bn_arcs$to[target_genes_list_indices[[i]]])
  }
  
  
  x <- data.frame()
  tf_tg_df <- data.frame()
  for(i in 1:length(target_genes_list_indices)){
    x <- bn_arcs[target_genes_list_indices[[i]],]
    tf_tg_df <- rbind(tf_tg_df,x)
  }
  
  
  
  #Get the details of which metabolites are getting affected by the Target genes
  
  #length(tf_tg_df$to) # 1049
  #length(unique(tf_tg_df$to)) # 783
  
  #load the gene subsytem information 
  #setwd("/data/users/cs20d300/WORK_MAIN/Integrated_ntwk/Ecoli/subsystems")
  #gene_subsys <- readRDS("Gene_subsystem_BIGGID_Names.rds")
  
  tf_tg_df[,ncol(tf_tg_df)+1] <- 0
  colnames(tf_tg_df)[ncol(tf_tg_df)] <- "tg_bigg_sym"
  
  for(i in 1:length(GS$Gene)){
    for(j in 1:length(tf_tg_df$to)){
      if(GS$gene_name[i]==tf_tg_df$to[j]){
        tf_tg_df$tg_bigg_sym[j]=GS$Gene[i]
      }
    }
  }
  
  #sum(tf_tg_df_xm201023$tg_bigg_sym==0) # number of 0s in the column of tg_bigg_sys - 374
  # 374 out of 1049 (length) was 0
  
  # subset of this modified tf_tg_df - to ensure genes that are not present in the metabolic model are removed  (2c-1c)
  # This is the dataframe to focus on for subsequent analysis (to calculate probabilities)
  tf_tg_df_met_model <- tf_tg_df[tf_tg_df$tg_bigg_sym !=0,]
  
  #length(tf_tg_df_met_model_xm201023$to)   # 675
  #length(unique(tf_tg_df_met_model_xm201023$to))   # 466
  
  
  tf_tg_df_met_model[,ncol(tf_tg_df_met_model)+1] <- 0
  colnames(tf_tg_df_met_model)[ncol(tf_tg_df_met_model)] <- "Regulation"
  
  for(i in 1:nrow(modified_interaction_df)){
    for(j in 1:nrow(tf_tg_df_met_model)){
      if(tf_tg_df_met_model$from[j]==modified_interaction_df$`TF-gene`[i]){
        tf_tg_df_met_model$Regulation[j]=modified_interaction_df$Interaction[i]
      }
    }
  }
  
  op_list <- list(tf_tg_df_met_model,modified_interaction_df)
  names(op_list) <- c("TF_TG_MET_df","TF_TG_int")
  
  
  
  #return(more_1_gene_int)
  #return(modified_interaction_df)
  return(op_list)}



Integration_part_2 <- function(op_lis_int_1, KO_gene, GS){
  
  get_the_evi_list <- function(op_lis){
    
    op_df_1 <- op_lis[[1]]
    op_df_2 <- op_lis[[2]]
    
    req_var <- unique(as.character(op_df_1[[1]]))
    extra_var <- op_df_2[[1]]
    
    needed_var <- intersect(req_var,extra_var)
    
    symb_vec <- list()
    for(j in 1:length(needed_var)){
      for(i in 1:nrow(op_df_2)){
        if(needed_var[j]==op_df_2[[1]][i]){
          if((op_df_2[[2]][i]=="-")==TRUE){
            symb_vec[j] <- paste("0")
          }else{
            symb_vec[j] <- paste("1")
          }
        }
      }
      
    }
    names(symb_vec) <- needed_var
    
    return(symb_vec)}
  
  
  e_l <- get_the_evi_list(op_lis_int_1)
  
  
  if(length(KO_gene)==0){
    e_l <- e_l
  }else{
    e_l <- append(e_l,"0")
    names(e_l)[length(e_l)] <- KO_gene
  }
  
  
  CP_cal_func <- function(op_lis, e){
    
    some_df_1 <- op_lis[[1]]
    some_df_2 <- op_lis[[2]]
    
    #e_l <- e
    
    tf_uni <- unique(some_df_1$from)
    tf_uni <- as.character(tf_uni)
    
    tg_uni <- unique(some_df_1$to)
    tg_uni <- as.character(tg_uni)
    
    new_TF_Int <- unique(some_df_2)
    
    tf_reg <- c()
    for(i in 1:length(tf_uni)){
      for(j in 1:nrow(new_TF_Int)){
        if(tf_uni[i]==new_TF_Int$`TF-gene`[j]){
          tf_reg <- c(tf_reg,new_TF_Int$Interaction[j])
        }
      }
    }
    
    cp_cal_tf_reg_df <- data.frame(tf_uni, tf_reg)  # table used as reference for cp evaluation 
    cp_cal_tf_reg_df[,ncol(cp_cal_tf_reg_df)+1] <- 0
    colnames(cp_cal_tf_reg_df)[ncol(cp_cal_tf_reg_df)] <- "tf_reg_num"
    
    for(i in 1:nrow(cp_cal_tf_reg_df)){
      if(cp_cal_tf_reg_df$tf_reg[i] == "+"){
        cp_cal_tf_reg_df$tf_reg_num[i] = 1
      }
    }
    
    
    reg_uni <- cp_cal_tf_reg_df$tf_reg_num
    
    
    all_tgs_eval_bde <- c()
    for(i in 1:length(tg_uni)){
      #set.seed(61)
      qtxt_1 <- paste0("cpquery(grn_PL, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
      x <- eval(parse(text = qtxt_1))
      all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
    }
    
    
    gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
    colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
    
    
    return(gene_cp_gg_df)}
  
  
  CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l)
  
  bid <- c() 
  
  for(i in 1:nrow(CP_DF_MET_MODEL)){
    for(j in 1:nrow(GS)){
      if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
        bid <- c(bid, GS$Gene[j])
      }
    }
  }
  
  CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
  colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
  
  int2_res <- list(CP_DF_MET_MODEL, e_l)
  names(int2_res) <- c("CP_Final","evidence_list")
  
  
  
  return(int2_res)}



############################################  step - 2d - evaluating upper bounds

to_get_upd_FVA <- function(gpr_vec, fva_r0_op){
  ### Getting the upper bounds of the reactions with modified GPR (probabilities)
  gpr_ids <- noquote(row.names(gpr_vec))
  gpr_ids <- as.numeric(gpr_ids)
  
  #FVA max value based on the GPR ids
  fva_max_ids <- fva_r0_op[[3]][gpr_ids]
  
  
  new_FVA_ub <- c()
  
  for(i in 1:nrow(gpr_vec)){
    mul <- fva_max_ids[i]*gpr_vec[[1]][i]
    new_FVA_ub <- c(new_FVA_ub, mul)
  }
  
  # new FVA table - update
  updated_FVA_table <- fva_r0_op
  updated_FVA_table$new_upper_bounds <- new_FVA_ub
  
  return(updated_FVA_table)}






########### To check for stoping criteria 

to_check_sink_rxn_iter <- function(fva_1, fva_2){
  
  sink_FVA_round_2 <- fva_2[[3]][6:8]
  
  round_2_vec <- c()
  for(i in 1:length(sink_FVA_round_2)){
    if(sink_FVA_round_2[i]>0){
      round_2_vec[i] = 1
    }else if(sink_FVA_round_2[i]==0){
      round_2_vec[i] = 0
    }
  }
  
  
  # Maximum fluxes for the sink reactions from round 2 of FVA
  sink_FVA_round_1 <- fva_1[[3]][6:8]
  
  round_1_vec <- c()
  for(i in 1:length(sink_FVA_round_1)){
    if(sink_FVA_round_1[i]>0){
      round_1_vec[i] = 1
    }else if(sink_FVA_round_1[i]==0){
      round_1_vec[i] = 0
    }
  }
  
  x_status <- identical(round_1_vec, round_2_vec)
  
  op <- list(x_status, round_1_vec, round_2_vec)
  names(op) <- c("x_status", "Round_1_FVA_max_bin", "Round_2_FVA_max_bin")
  
  return(op)}





#######################################################################################################################
########################  Application for the same 
library(readxl) 
library(matlabr)
library(bnlearn)

#### Initialization step 

curr_wd <- c("D:/work/Integrated_network_model/Toy_model/auto_res")


## GRN - module
setwd("D:/work/Integrated_network_model/Toy_model/results/GRN_module/")
grn_SL <- readRDS("structure_learning.RDS")
grn_PL <- readRDS("parameter_learning.RDS")


## Metabolic module

setwd(curr_wd)
matlabr::run_matlab_script("MAT_1b_Initialization.m", display = TRUE, verbose = TRUE)


FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
FVA_round_0 <- as.data.frame(FVA_round_0)


#### Iteration step 

## get all the variables loaded
setwd("D:/work/Integrated_network_model/Toy_model/data/")
met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
gene_subsys <- readRDS("gene_subsys.RDS")


FVA_up <- FVA_round_0

count <- 1

repeat{
  
  FVA_to_check <- readxl::read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
  colnames(FVA_to_check) <- c("Reactions", "Minimum_flux", "Maximum_flux")
  FVA_to_check <- as.data.frame(FVA_to_check)
  
  FVA_i <- FVA_to_check
  
  # step - 2b
  new_met_gene_reg_data <- get_the_binary_data(met_gene_reg_data,FVA_i)
  # step - 2c - part 1
  op_intg_1 <- Integration_part_1(new_met_gene_reg_data, grn_SL, gene_subsys)
  # step - 2c - part 2
  op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = c(), gene_subsys)
  #write the o/p from 2c to xlsx for matlab function 2
  writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
  
  # step - 2d - part 1
  setwd(curr_wd)
  matlabr::run_matlab_script("MAT_2d_part_1_Iteration.m", display = TRUE, verbose = TRUE)
  
  gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
  gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
  colnames(gpr_eval_round_i) <- "gpr_eval"
  
  # step - 2d - part 2
  new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
  
  writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
  
  # step - 2a - run FVA with 0 objective
  setwd(curr_wd)
  matlabr::run_matlab_script("MAT_2a_FVA_Iteration.m", display = TRUE, verbose = TRUE)
  
  
  FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
  FVA_iplus1 <- as.data.frame( FVA_iplus1)
  colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
  
  
  iter_op <- to_check_sink_rxn_iter(FVA_i, FVA_iplus1)
  count <- count + 1
  
  # step - 2e - break
  if(iter_op[[1]]==TRUE){
    break
  }
  
}





























