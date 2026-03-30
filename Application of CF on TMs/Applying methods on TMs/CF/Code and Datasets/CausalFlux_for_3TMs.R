##################################################################################
curr_wd <- c("D:/work/Integrated_network_model/Git_hub_codes/TMS/")


####################################################### TM codes

TM1_Single_KO_CF_S <- function(curr_wd,pe,mi, xun){
  
  GRN_CS <- function(KO,sl, pl,g, op1){
    
    if( 0 %in% length(op1) & "WT" %in% KO){
      yu <- list(sl,pl)
      names(yu) <- c("CS_SL", "CS_PL")
    }else{
      
      if(length(op1) != 0){
        to_d <- op1[[2]]
        
        bh <- vector(mode = "list", length = nrow(to_d))
        names(bh) <- to_d[[1]]
        
        for(i in 1:length(bh)){
          bh[[i]] <- to_d[[3]][i]
        }
        
        if(KO == "WT"){
          bh <- bh
        }else{
          nx <- names(bh)
          if(KO %in% nx){
            
            koi <- which(nx %in% KO)
            bh[[koi]] <- 0
            
          }else{
            bh[[length(bh)+1]] <- 0
            names(bh)[length(bh)] <- KO
          }
          
        }
        
      }else{
        bh <- list()
        bh[[1]] <- 0
        names(bh) <- KO
      }
      
      
      
      #bo <- paste0("mutilated(sl, evidence = list(",noquote(paste(KO)),"= 0))")
      
      #bo_s <- eval(parse(text = bo))
      bo_s <- mutilated(sl,evidence = bh)
      
      bo_p <- bn.fit(bo_s, g, method = "bayes")
      
      if(length(op1) != 0){
        bess <- c(to_d[[1]], KO)
        
        keep_ind <- which(names(bo_p) %in% bess)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        
        
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        
        cptB = matrix(c(0, 1), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_2 <- as.table(cptB)
        
        
        # Replace contents of A with B except at the chosen index
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        
        for(i in 1:length(keep_ind_order)){
          xox <- which(names(bh) %in% keep_ind_order[i])
          
          if(bh[[xox]]==1){
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptB)", sep = "")  
            
            eval(parse(text = so))
          }else{
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptA)", sep = "")  
            
            eval(parse(text = so))
          }
          
        }
        
        
        # bo_p_modified[-keep_ind] <- pl[-keep_ind]
        # 
        # # Modify the content at the unchanged index
        # #bo_p_modified$A <- as.table(cptA)
        # 
        # so <- paste("bo_p_modified$",KO," <- as.table(cptA)", sep = "")  
        # 
        # eval(parse(text = so))
        # 
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
      }else{
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        keep_ind <- which(names(bo_p) %in% KO)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        so <- paste("bo_p_modified$",keep_ind_order," <- as.table(cptA)", sep = "")  
        
        eval(parse(text = so))
        
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
        
        
      }
      
    }
    
    
    
    
    
    return(yu)}
  
  
  get_the_binary_data <- function(df1, df2, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>0){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    
    grnarcs <- data.frame(sl$arcs)
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    x <- setdiff(nt[[3]], ft)
    
    xid <- which(nt[[3]] %in% x)
    
    df1 <- nt[-xid,]
    
    
    return(df1)}
  
  
  # df1 <- met_gene_reg_data
  # df2 <- FVA_i
  # df3 <- FVA_XP
  # p <- percen
  # sl <- grn_SL
  
  get_the_binary_data_iter <- function(df1, df2, df3, p, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>p*df3[[3]][j]){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    grnarcs <- data.frame(sl$arcs)
    
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    #x <- setdiff(nt[[3]], ft)
    #x <- setdiff(ft,nt[[3]])
    #xid <- which(nt[[3]] %in% x)
    xid <- which(nt[[3]] %in% ft)
    
    df1 <- nt[xid,]
    
    
    return(df1)}
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
  # df_1 <- NEW_MGR
  # bn_obj <- grn_SL
  # GS <- gene_subsys
  # KO <- xu_n[j]
  
  Integration_part_1 <- function(df_1, bn_obj, GS, KO){
    all <- df_1[[6]]
    
    if((all(all == 0))== FALSE){
      # removing the unkown interactions
      #new_met_gene_reg_data #df1
      q_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[4]][i]=="?"
        if("?" %in% df_1[[4]][i]){
          q_idx <- c(q_idx,i)
        }
      }
      
      p_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[6]][i]=="0"
        if("0" %in% df_1[[6]][i]){
          p_idx <- c(p_idx,i)
        }
      }
      
      pq_idx <- union(q_idx, p_idx)
      #df_1 <- new_met_gene_reg_data
      if(length(pq_idx)!=0){
        df_1 <- df_1[-pq_idx,]
      }else{
        df_1 <- df_1
      }
      
      
      # to remove those TFs with both "+" and "-" interactions
      # get the unique TFs
      unique_tf <- unique(df_1[[3]])
      
      
      if(length(unique_tf)==length(df_1$`TF-gene`)){
        
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
        
        # more_1_gene <- c()
        # for(i in 1:length(x_id)){
        #   more_1_gene <- c(more_1_gene,names(UNI_TF_name_ids[x_id[i]]))
        # }
        
        # more_1_gene_int <- list()
        # for(i in 1:length(x_id)){
        #   more_1_gene_int[[i]] <- df_1[[4]][UNI_TF_name_ids[[x_id[i]]]]
        # }
        
        # names(more_1_gene_int) <- more_1_gene
        
        # to get a df with only 1 rep of regulation
        only_1_idx <- c()
        for(i in 1:length(UNI_TF_name_ids)){
          if(length(UNI_TF_name_ids[[i]])==1){
            only_1_idx <- c(only_1_idx,i)
          }
        }
        
        req_names <- names(UNI_TF_name_ids)[only_1_idx]
        
        modified_interaction_df <- df_1[which(df_1[[3]] %in% req_names),3:4]
        
        
      } else if (length(unique_tf)!=length(df_1$`TF-gene`)){
        
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
          if("0" %in% df_1[[6]][i]){
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
        
        
      }
      
      
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
      tf_tg_df[,ncol(tf_tg_df)+1] <- 0
      colnames(tf_tg_df)[ncol(tf_tg_df)] <- "tg_bigg_sym"
      
      for(i in 1:length(GS$Gene)){
        for(j in 1:length(tf_tg_df$to)){
          if(GS$gene_name[i]==tf_tg_df$to[j]){
            tf_tg_df$tg_bigg_sym[j]=GS$Gene[i]
          }
        }
      }
      
      
      # subset of this modified tf_tg_df - to ensure genes that are not present in the metabolic model are removed  (2c-1c)
      # This is the dataframe to focus on for subsequent analysis (to calculate probabilities)
      tf_tg_df_met_model <- tf_tg_df[tf_tg_df$tg_bigg_sym !=0,]
      
      
      tf_tg_df_met_model[,ncol(tf_tg_df_met_model)+1] <- 0
      colnames(tf_tg_df_met_model)[ncol(tf_tg_df_met_model)] <- "Regulation"
      
      for(i in 1:nrow(modified_interaction_df)){
        for(j in 1:nrow(tf_tg_df_met_model)){
          if(tf_tg_df_met_model$from[j]==modified_interaction_df$`TF-gene`[i]){
            tf_tg_df_met_model$Regulation[j]=modified_interaction_df$Interaction[i]
          }
        }
      }
      
      
      
      
      biis <- c()
      for(i in 1:nrow(modified_interaction_df)){
        if(modified_interaction_df$Interaction[i] == "+"){
          biis[i] <- 1
        }else{
          biis[i] <- 0
        }
      }
      modified_interaction_df$status <- biis
      
      y <- union(bn_arcs$from, bn_arcs$to)
      
      y_ind <- which(modified_interaction_df$`TF-gene` %in% y)
      
      modified_interaction_df <- modified_interaction_df[y_ind,]
      
      
      if(KO == "WT"){
        final_df <- tf_tg_df_met_model
      }else{
        bn_arcs <- bn_obj[[3]]
        bn_arcs <- as.data.frame(bn_arcs)
        
        bn_f <- which(bn_arcs$from %in% KO)
        
        
        if(length(bn_f) != 0){
          
          bn_f_a <- bn_arcs[bn_f,]
          bn_f_a[,ncol(bn_f_a)+1] <- 0
          colnames(bn_f_a)[ncol(bn_f_a)] <- "tg_bigg_sym"
          
          for(i in 1:length(GS$Gene)){
            for(j in 1:length(bn_f_a$to)){
              if(GS$gene_name[i]==bn_f_a$to[j]){
                bn_f_a$tg_bigg_sym[j]=GS$Gene[i]
              }
            }
          }
          
          bn_f_a$Regulation <- rep("-", nrow(bn_f_a))
          
          x0 <- which(bn_f_a[[3]] %in% 0)
          bn_f_a <- bn_f_a[-x0,]
          
          
          final_df <- rbind(tf_tg_df_met_model, bn_f_a)
          
          final_df <- unique(final_df)
        }else{
          final_df <- tf_tg_df_met_model
        }
        
      }
      
      
      op_list <- list(final_df,modified_interaction_df)
      names(op_list) <- c("TF_TG_MET_df","TF_TG_int")
      
    }else{
      op_list <- c()
    }
    
    
    
    return(op_list)}
  
  Integration_initial <- function(KO_gene, GS,stra,para, g ,wds, ct,fu){
    if("WT" %in% KO_gene){
      setwd(wds)
      matlabr::run_matlab_script("MAT_1b_Initialization_TOY.m", display = TRUE, verbose = TRUE)
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      # 
      
      # KO_gene = "gnd"
      # stra = grn_SL
      # para = grn_PL
      # g <- gge
      # 
      
      
    }else{
      
      x <- GRN_CS(KO_gene,stra,para, g,c())
      xcspl <- x$CS_PL
      
      ga <- data.frame(stra$arcs)
      
      tg_g <- c()
      x_l <- which(ga$from %in% KO_gene)
      if(length(x_l) != 0){
        tg_g <- ga$to[x_l]
        tg_g[length(tg_g)+1] <- KO_gene
      }else{
        tg_g <- KO_gene
      }
      
      
      all_tgs_eval_bde <- c()
      #set.seed(61)
      for(i in 1:length(tg_g)){
        ##   should b before ( and after e
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
        qtxt_1 <- paste0("cpquery(xcspl, ","(",paste(noquote(tg_g[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
        x <- eval(parse(text = qtxt_1))
        all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
      }
      
      
      gene_cp_gg_df <- data.frame(tg_g,all_tgs_eval_bde)
      colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if(gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      # colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      # 
      bid <- c() 
      
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if((gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]) == TRUE){
      #       #bid <- c(bid, GS$Gene[j])
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(gene_cp_gg_df$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      
      gene_cp_gg_df[is.na(gene_cp_gg_df)] <- 0
      
      int2_res <- list(gene_cp_gg_df)
      names(int2_res) <- c("CP_Final")
      
      
      writexl::write_xlsx(int2_res$CP_Final, paste0(wds,"/CP_round_i.xlsx"),col_names = TRUE)
      
      cp <- read_xlsx(paste0(wds,"/CP_round_i.xlsx"), col_names = FALSE)
      write.csv(cp, paste0("CP_incorp_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY.m", display = TRUE, verbose = TRUE)
      
      gpr_eval_round_i <- read_xlsx(paste0(wds,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
      gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
      colnames(gpr_eval_round_i) <- "gpr_eval"
      
      # step - 2d - part 2
      new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, fu)
      
      writexl::write_xlsx(new_upd_fva,paste0(wds,"/Updated_FVA_round_i.xlsx"))
      
      mup <- read_xlsx(paste0(wds,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
      write.csv(mup, paste0("Updated_FVA_round_incorp_",ct,".csv"))
      
      
      # step - 2a - run FVA with 0 objective
      setwd(wds)
      matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY.m", display = TRUE, verbose = TRUE)
      
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      
    }
    return(FVA_incorp)}
  
  #Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
  
  
  Integration_part_2 <- function(op_lis_int_1, KO_gene, GS,para, strat){
    
    if(KO_gene == "WT" & length(op_lis_int_1)==0){
      vm <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
      colnames(vm) <- c("MET_MODEL_TG_genes", "Probability", "Bigg_symb")
      int2_res <- vector("list",2)
      int2_res[[1]] <- vm
      names(int2_res) <- c("CP_Final","evidence_list")
    }else if(KO_gene == "WT"){
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        some_df_1 <- op_lis[[1]]
        some_df_2 <- op_lis[[2]]
        
        ga <- data.frame(stra$arcs)
        
        # tf_uni <- unique(some_df_1$from)
        # tf_uni <- as.character(tf_uni)
        
        tg_uni <- unique(some_df_1$to)
        tg_uni <- as.character(tg_uni)
        
        
        tf_uni  <- some_df_1$from
        tf_uni <- as.character(tf_uni)
        
        sta_vec <- c()
        for(i in 1:nrow(some_df_1)){
          if(some_df_1$Regulation[i]== "+"){
            sta_vec[i] <- 1
          }else{
            sta_vec[i] <- 0
          }
        }
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(tf_uni[i])),"==",paste(noquote(sta_vec[i])),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l, stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
      
    }else{
      
      # get_the_evi_list <- function(op_lis){
      #   
      #   op_df_1 <- op_lis[[1]]
      #   op_df_2 <- op_lis[[2]]
      #   
      #   req_var <- unique(as.character(op_df_1[[1]]))
      #   extra_var <- op_df_2[[1]]
      #   
      #   needed_var <- intersect(req_var,extra_var)
      #   
      #   symb_vec <- list()
      #   for(j in 1:length(needed_var)){
      #     symb_vec[j] <- paste("1")
      #     
      #   }
      #   names(symb_vec) <- needed_var
      #   
      #   return(symb_vec)}
      # 
      # 
      # e_l <- get_the_evi_list(op_lis_int_1)
      # 
      # 
      # if(length(KO_gene)==0){
      #   e_l <- e_l
      # } else if(KO_gene %in% names(e_l)){
      #   x <- which(names(e_l) %in% KO_gene)
      #   e_l[[x]] <- "0"
      # }else {
      #   e_l <- append(e_l,"0")
      #   names(e_l)[length(e_l)] <- KO_gene
      # }
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        if(length(op_lis) != 0){
          some_df_1 <- op_lis[[1]]
          some_df_2 <- op_lis[[2]]
          
          #e_l <- e
          
          tf_uni <- unique(some_df_1$from)
          tf_uni <- as.character(tf_uni)
          
          tg_uni <- unique(some_df_1$to)
          tg_uni <- as.character(tg_uni)
          tg_uni[length(tg_uni)+1] <- KO_gene
          
        }else{
          ga <- data.frame(stra$arcs)
          
          tg_uni <- c()
          x_l <- which(ga$from %in% KO_gene)
          if(length(x_l) != 0){
            tg_uni <- ga$to[x_l]
          }else{
            tg_uni <- KO_gene
          }
        }
        
        
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l,stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
    }
    
    
    
    
    
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
    
    
    ############################### Lower bounds modification part
    gv <- gpr_vec[[1]]
    xo <- fva_r0_op[[2]]
    
    xo_ind <- which(xo < 0)
    
    for(i in 1:length(xo_ind)){
      xo[xo_ind[i]] <-  xo[xo_ind[i]]*gv[xo_ind[i]]
    }
    
    
    
    # new FVA table - update
    updated_FVA_table <- fva_r0_op
    updated_FVA_table$new_upper_bounds <- new_FVA_ub
    updated_FVA_table$new_lower_bounds <- xo
    
    return(updated_FVA_table)}
  
  
  ########### To check for stoping criteria 
  
  to_check_sink_rxn_iter <- function(fva_pre, fva_1, fva_2,r_begin, r_end,p){
    
    sink_FVA_round_2 <- fva_2[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    
    sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] 
    
    sink_FVA_round_pre <- fva_pre[[3]][r_begin:r_end] 
    
    round_2_vec <- c()
    
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }
    
    round_1_vec <- c()
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }
    
    
    
    # 
    # for(i in 1:length(sink_FVA_round_2)){
    #   if(sink_FVA_round_2[i]> p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 1
    #   }else if(sink_FVA_round_2[i]< p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 0
    #   }
    # }
    
    
    # Maximum fluxes for the sink reactions from round 2 of FVA
    
    # sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    # 
    # round_1_vec <- c()
    # for(i in 1:length(sink_FVA_round_1)){
    #   if(sink_FVA_round_1[i]>0){
    #     round_1_vec[i] = 1
    #   }else if(sink_FVA_round_1[i]<0){
    #     round_1_vec[i] = 0
    #   }
    # }
    
    x_status <- identical(round_1_vec, round_2_vec)
    
    op <- list(x_status, round_1_vec, round_2_vec)
    names(op) <- c("x_status", "Round_1_FVA_max_bin", "Round_2_FVA_max_bin")
    
    return(op)}
  
  
  #######################################################################################################################
  # ########################  Application for the same 
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  library(writexl)
  library(tidyverse)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    
    setwd(curr_wd)
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_3200/"))
    }
    
    
    grn_SL <- readRDS("Structure_learning.rds")
    
    saveRDS(grn_SL, file = paste0("grn_sl_",count-1,".RDS"))
    
    grn_PL <- readRDS("Parameter_learning.rds")
    
    saveRDS(grn_PL, file = paste0("grn_pl_",count-1,".RDS"))
    
    gge <- read.csv(paste0("Bin_GE_TM1_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    
    setwd(curr_wd)
    matlabr::run_matlab_script("MAT_1b_Initialization_TOY.m", display = TRUE, verbose = TRUE)
    
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM1/grn_3200/"))
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    count <- 1
    
    setwd(curr_wd)
    FVA_to_check <- readxl::read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    colnames(FVA_to_check) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_to_check <- as.data.frame(FVA_to_check)
    
    write.csv(FVA_to_check, paste0("FVA_incorp_",count-1,".csv"))
    
    FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_incorp, paste0("FBA_incorp_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_to_check
    
    FVA_III <- FVA_to_check
    #FVA_III <- FBA_incorp
    
    FVA_II <- FVA_III
    
    xu_n <- xun
    
    
    for(j in 1:length(xu_n)){
      
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_to_check
      
      count <- 1
      
      
      FVA_III <- FVA_to_check
      #FVA_III <- FBA_incorp
      
      FVA_II <- FVA_III
      
      
      FVA_one <- FVA_II
      
      ifva <- Integration_initial(xu_n[j], gene_subsys, grn_SL, grn_PL,gge, curr_wd,count, FVA_up)
      
      
      FVA_P <- ifva 
      #FVA_P <- FBA_incorp
      
      FVA_prev <- FVA_III 
      
      repeat{
        
        #NEW_MGR <- NEW_MGR_prime
        FVA_XP <- FVA_prev
        
        FVA_q <- FVA_P
        
        FVA_i <- FVA_q
        
        NEW_MGR <- get_the_binary_data_iter(met_gene_reg_data,FVA_i,FVA_XP, percen, grn_SL)
        
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(NEW_MGR, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        write.csv(op_intg_1[[1]], paste0("TF_TG_MET_",count+1,".csv"))
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_incorp_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        mup <- read_xlsx(paste0(curr_wd,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
        write.csv(mup, paste0("Updated_FVA_round_incorp_",count+1,".csv"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY.m", display = TRUE, verbose = TRUE)
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        FVA_iplus1 <- as.data.frame( FVA_iplus1)
        colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        
        #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_incorp_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_incorp_",count+1,".csv"))
        
        
        iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FVA_iplus1,6,8, percen)
        #iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FBA_iplus1,2713,2778, percen)
        
        
        write.csv(iter_op$Round_1_FVA_max_bin,paste0("BV_",count,".csv"))
        write.csv(iter_op$Round_2_FVA_max_bin,paste0("BV_",count+1,".csv"))
        
        
        
        count <- count + 1
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        #NEW_MGR_prime <- iter_op$Round_2_FVA_max_bin
        
        FVA_prev <- FVA_i
        
        FVA_P <- FVA_iplus1
        
      }
      
      bv <- c()
      for(i in 0:(count)){
        bv <- c(bv,paste0("BV_",i,".csv"))
      }
      
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_incorp_",i,".csv"))
      }
      
      cc1 <- c()
      for(i in 0:(count)){
        cc1 <- c(cc1,paste0("grn_sl_",i,".RDS"))
      }
      
      cc2 <- c()
      for(i in 0:(count)){
        cc2 <- c(cc2,paste0("grn_pl_",i,".RDS"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_incorp_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_incorp_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      mm <- c()
      for(i in 1:(count)){
        mm <- c(mm,paste0("TF_TG_MET_",i,".csv"))
      }
      
      my_files <- c("FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",vv,bb,pp,uu,mm,cc1,cc2,bv)
      curr_wd_ <- paste0(curr_wd,"/")
      ko_dir <- paste0(curr_wd_,xu_n[j],"/")
      
      dir.create(paste0(curr_wd_,xu_n[j]))
      
      file.rename(from = paste0(curr_wd_, my_files),to = paste0(ko_dir, my_files))
      
    }
    library(ff)
    dir.create(paste0(curr_wd,"/KO_data_TM1_",ee[e]))
    from <- curr_wd            #Current path of your folder
    to   <- curr_wd            #Path you want to move it.
    
    m <- xu_n
    
    for(i in 1:length(m)){
      path1 <- paste0(from,"/",m[i])
      path2 <- paste0(to,"/KO_data_TM1_",ee[e],"/",m[i])
      file.rename(path1,path2)
      
    }
  }
  
}




TM2_Single_KO_CF_S <- function(curr_wd,pe,mi, xun){
  
  # KO <- xu_n[j]
  # sl <- grn_SL
  # pl <- grn_PL
  # g <- gge
  # op1 <- op_intg_1
  
  GRN_CS <- function(KO,sl, pl,g, op1){
    
    if(length(op1) == 0 & KO == "WT"){
      yu <- list(sl,pl)
      names(yu) <- c("CS_SL", "CS_PL")
    }else{
      
      if(length(op1) != 0){
        to_d <- op1[[2]]
        
        bh <- vector(mode = "list", length = nrow(to_d))
        names(bh) <- to_d[[1]]
        
        for(i in 1:length(bh)){
          bh[[i]] <- to_d[[3]][i]
        }
        
        if(KO == "WT"){
          bh <- bh
        }else{
          nx <- names(bh)
          if(KO %in% nx){
            
            koi <- which(nx %in% KO)
            bh[[koi]] <- 0
            
          }else{
            bh[[length(bh)+1]] <- 0
            names(bh)[length(bh)] <- KO
          }
          
        }
        
      }else{
        bh <- list()
        bh[[1]] <- 0
        names(bh) <- KO
      }
      
      
      
      #bo <- paste0("mutilated(sl, evidence = list(",noquote(paste(KO)),"= 0))")
      
      #bo_s <- eval(parse(text = bo))
      bo_s <- mutilated(sl,evidence = bh)
      
      bo_p <- bn.fit(bo_s, g, method = "bayes")
      
      if(length(op1) != 0){
        bess <- c(to_d[[1]], KO)
        
        keep_ind <- which(names(bo_p) %in% bess)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        
        
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        
        cptB = matrix(c(0, 1), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_2 <- as.table(cptB)
        
        
        # Replace contents of A with B except at the chosen index
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        
        for(i in 1:length(keep_ind_order)){
          xox <- which(names(bh) %in% keep_ind_order[i])
          
          if(bh[[xox]]==1){
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptB)", sep = "")  
            
            eval(parse(text = so))
          }else{
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptA)", sep = "")  
            
            eval(parse(text = so))
          }
          
        }
        
        
        # bo_p_modified[-keep_ind] <- pl[-keep_ind]
        # 
        # # Modify the content at the unchanged index
        # #bo_p_modified$A <- as.table(cptA)
        # 
        # so <- paste("bo_p_modified$",KO," <- as.table(cptA)", sep = "")  
        # 
        # eval(parse(text = so))
        # 
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
      }else{
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        keep_ind <- which(names(bo_p) %in% KO)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        so <- paste("bo_p_modified$",keep_ind_order," <- as.table(cptA)", sep = "")  
        
        eval(parse(text = so))
        
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
        
        
      }
      
    }
    
    
    
    
    
    return(yu)}
  
  get_the_binary_data <- function(df1, df2, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>0){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    
    grnarcs <- data.frame(sl$arcs)
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    x <- setdiff(nt[[3]], ft)
    
    xid <- which(nt[[3]] %in% x)
    
    df1 <- nt[-xid,]
    
    
    return(df1)}
  
  
  # df1 <- met_gene_reg_data
  # df2 <- FVA_i
  # df3 <- FVA_XP
  # p <- percen
  # sl <- grn_SL
  
  get_the_binary_data_iter <- function(df1, df2, df3, p, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>p*df3[[3]][j]){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    grnarcs <- data.frame(sl$arcs)
    
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    #x <- setdiff(nt[[3]], ft)
    #x <- setdiff(ft,nt[[3]])
    #xid <- which(nt[[3]] %in% x)
    xid <- which(nt[[3]] %in% ft)
    
    df1 <- nt[xid,]
    
    
    return(df1)}
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
  # df_1 <- NEW_MGR
  # bn_obj <- grn_SL
  # GS <- gene_subsys
  # KO <- xu_n[j]
  
  Integration_part_1 <- function(df_1, bn_obj, GS, KO){
    all <- df_1[[6]]
    
    if((all(all == 0))== FALSE){
      # removing the unkown interactions
      #new_met_gene_reg_data #df1
      q_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[4]][i]=="?"
        if("?" %in% df_1[[4]][i]){
          q_idx <- c(q_idx,i)
        }
      }
      
      p_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[6]][i]=="0"
        if("0" %in% df_1[[6]][i]){
          p_idx <- c(p_idx,i)
        }
      }
      
      pq_idx <- union(q_idx, p_idx)
      #df_1 <- new_met_gene_reg_data
      if(length(pq_idx)!=0){
        df_1 <- df_1[-pq_idx,]
      }else{
        df_1 <- df_1
      }
      
      
      # to remove those TFs with both "+" and "-" interactions
      # get the unique TFs
      unique_tf <- unique(df_1[[3]])
      
      
      if(length(unique_tf)==length(df_1$`TF-gene`)){
        
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
        
        # more_1_gene <- c()
        # for(i in 1:length(x_id)){
        #   more_1_gene <- c(more_1_gene,names(UNI_TF_name_ids[x_id[i]]))
        # }
        
        # more_1_gene_int <- list()
        # for(i in 1:length(x_id)){
        #   more_1_gene_int[[i]] <- df_1[[4]][UNI_TF_name_ids[[x_id[i]]]]
        # }
        
        # names(more_1_gene_int) <- more_1_gene
        
        # to get a df with only 1 rep of regulation
        only_1_idx <- c()
        for(i in 1:length(UNI_TF_name_ids)){
          if(length(UNI_TF_name_ids[[i]])==1){
            only_1_idx <- c(only_1_idx,i)
          }
        }
        
        req_names <- names(UNI_TF_name_ids)[only_1_idx]
        
        modified_interaction_df <- df_1[which(df_1[[3]] %in% req_names),3:4]
        
        
      } else if (length(unique_tf)!=length(df_1$`TF-gene`)){
        
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
        
        
      }
      
      
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
      tf_tg_df[,ncol(tf_tg_df)+1] <- 0
      colnames(tf_tg_df)[ncol(tf_tg_df)] <- "tg_bigg_sym"
      
      for(i in 1:length(GS$Gene)){
        for(j in 1:length(tf_tg_df$to)){
          if(GS$gene_name[i]==tf_tg_df$to[j]){
            tf_tg_df$tg_bigg_sym[j]=GS$Gene[i]
          }
        }
      }
      
      
      # subset of this modified tf_tg_df - to ensure genes that are not present in the metabolic model are removed  (2c-1c)
      # This is the dataframe to focus on for subsequent analysis (to calculate probabilities)
      tf_tg_df_met_model <- tf_tg_df[tf_tg_df$tg_bigg_sym !=0,]
      
      
      tf_tg_df_met_model[,ncol(tf_tg_df_met_model)+1] <- 0
      colnames(tf_tg_df_met_model)[ncol(tf_tg_df_met_model)] <- "Regulation"
      
      for(i in 1:nrow(modified_interaction_df)){
        for(j in 1:nrow(tf_tg_df_met_model)){
          if(tf_tg_df_met_model$from[j]==modified_interaction_df$`TF-gene`[i]){
            tf_tg_df_met_model$Regulation[j]=modified_interaction_df$Interaction[i]
          }
        }
      }
      
      
      
      
      biis <- c()
      for(i in 1:nrow(modified_interaction_df)){
        if(modified_interaction_df$Interaction[i] == "+"){
          biis[i] <- 1
        }else{
          biis[i] <- 0
        }
      }
      modified_interaction_df$status <- biis
      
      y <- union(bn_arcs$from, bn_arcs$to)
      
      y_ind <- which(modified_interaction_df$`TF-gene` %in% y)
      
      modified_interaction_df <- modified_interaction_df[y_ind,]
      
      
      if(KO == "WT"){
        final_df <- tf_tg_df_met_model
      }else{
        bn_arcs <- bn_obj[[3]]
        bn_arcs <- as.data.frame(bn_arcs)
        
        bn_f <- which(bn_arcs$from %in% KO)
        
        
        if(length(bn_f) != 0){
          
          bn_f_a <- bn_arcs[bn_f,]
          bn_f_a[,ncol(bn_f_a)+1] <- 0
          colnames(bn_f_a)[ncol(bn_f_a)] <- "tg_bigg_sym"
          
          for(i in 1:length(GS$Gene)){
            for(j in 1:length(bn_f_a$to)){
              if(GS$gene_name[i]==bn_f_a$to[j]){
                bn_f_a$tg_bigg_sym[j]=GS$Gene[i]
              }
            }
          }
          
          bn_f_a$Regulation <- rep("-", nrow(bn_f_a))
          
          x0 <- which(bn_f_a[[3]] %in% 0)
          bn_f_a <- bn_f_a[-x0,]
          
          
          final_df <- rbind(tf_tg_df_met_model, bn_f_a)
          
          final_df <- unique(final_df)
        }else{
          final_df <- tf_tg_df_met_model
        }
        
      }
      
      
      op_list <- list(final_df,modified_interaction_df)
      names(op_list) <- c("TF_TG_MET_df","TF_TG_int")
      
    }else{
      op_list <- c()
    }
    
    
    
    return(op_list)}
  
  Integration_initial <- function(KO_gene, GS,stra,para, g ,wds, ct,fu){
    if(KO_gene == "WT"){
      setwd(wds)
      matlabr::run_matlab_script("MAT_1b_Initialization_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      # 
      
      # KO_gene = "gnd"
      # stra = grn_SL
      # para = grn_PL
      # g <- gge
      # 
      
      
    }else{
      
      x <- GRN_CS(KO_gene,stra,para, g,c())
      xcspl <- x$CS_PL
      
      saveRDS(x$CS_SL, file = paste0("grn_sl_",ct,".RDS"))
      saveRDS(x$CS_PL, file = paste0("grn_pl_",ct,".RDS"))
      
      ga <- data.frame(stra$arcs)
      
      tg_g <- c()
      x_l <- which(ga$from %in% KO_gene)
      if(length(x_l) != 0){
        tg_g <- ga$to[x_l]
        tg_g[length(tg_g)+1] <- KO_gene
      }else{
        tg_g <- KO_gene
      }
      
      
      all_tgs_eval_bde <- c()
      #set.seed(61)
      for(i in 1:length(tg_g)){
        ##   should b before ( and after e
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
        qtxt_1 <- paste0("cpquery(xcspl, ","(",paste(noquote(tg_g[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
        x <- eval(parse(text = qtxt_1))
        all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
      }
      
      
      gene_cp_gg_df <- data.frame(tg_g,all_tgs_eval_bde)
      colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if(gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      # colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      # 
      bid <- c() 
      
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if((gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]) == TRUE){
      #       #bid <- c(bid, GS$Gene[j])
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(gene_cp_gg_df$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      
      gene_cp_gg_df[is.na(gene_cp_gg_df)] <- 0
      
      int2_res <- list(gene_cp_gg_df)
      names(int2_res) <- c("CP_Final")
      
      
      writexl::write_xlsx(int2_res$CP_Final, paste0(wds,"/CP_round_i.xlsx"),col_names = TRUE)
      
      cp <- read_xlsx(paste0(wds,"/CP_round_i.xlsx"), col_names = FALSE)
      write.csv(cp, paste0("CP_incorp_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
      
      gpr_eval_round_i <- read_xlsx(paste0(wds,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
      gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
      colnames(gpr_eval_round_i) <- "gpr_eval"
      
      # step - 2d - part 2
      new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, fu)
      
      writexl::write_xlsx(new_upd_fva,paste0(wds,"/Updated_FVA_round_i.xlsx"))
      
      mup <- read_xlsx(paste0(wds,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
      write.csv(mup, paste0("Updated_FVA_round_incorp_",ct,".csv"))
      
      
      # step - 2a - run FVA with 0 objective
      setwd(wds)
      matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
      
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      
    }
    return(FVA_incorp)}
  
  #Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
  
  Integration_part_2 <- function(op_lis_int_1, KO_gene, GS,para, strat){
    
    if(KO_gene == "WT" & length(op_lis_int_1)==0){
      vm <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
      colnames(vm) <- c("MET_MODEL_TG_genes", "Probability", "Bigg_symb")
      int2_res <- vector("list",2)
      int2_res[[1]] <- vm
      names(int2_res) <- c("CP_Final","evidence_list")
    }else if(KO_gene == "WT"){
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        some_df_1 <- op_lis[[1]]
        some_df_2 <- op_lis[[2]]
        
        ga <- data.frame(stra$arcs)
        
        # tf_uni <- unique(some_df_1$from)
        # tf_uni <- as.character(tf_uni)
        
        tg_uni <- unique(some_df_1$to)
        tg_uni <- as.character(tg_uni)
        
        
        tf_uni  <- some_df_1$from
        tf_uni <- as.character(tf_uni)
        
        sta_vec <- c()
        for(i in 1:nrow(some_df_1)){
          if(some_df_1$Regulation[i]== "+"){
            sta_vec[i] <- 1
          }else{
            sta_vec[i] <- 0
          }
        }
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(tf_uni[i])),"==",paste(noquote(sta_vec[i])),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l, stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
      
    }else{
      
      # get_the_evi_list <- function(op_lis){
      #   
      #   op_df_1 <- op_lis[[1]]
      #   op_df_2 <- op_lis[[2]]
      #   
      #   req_var <- unique(as.character(op_df_1[[1]]))
      #   extra_var <- op_df_2[[1]]
      #   
      #   needed_var <- intersect(req_var,extra_var)
      #   
      #   symb_vec <- list()
      #   for(j in 1:length(needed_var)){
      #     symb_vec[j] <- paste("1")
      #     
      #   }
      #   names(symb_vec) <- needed_var
      #   
      #   return(symb_vec)}
      # 
      # 
      # e_l <- get_the_evi_list(op_lis_int_1)
      # 
      # 
      # if(length(KO_gene)==0){
      #   e_l <- e_l
      # } else if(KO_gene %in% names(e_l)){
      #   x <- which(names(e_l) %in% KO_gene)
      #   e_l[[x]] <- "0"
      # }else {
      #   e_l <- append(e_l,"0")
      #   names(e_l)[length(e_l)] <- KO_gene
      # }
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        if(length(op_lis) != 0){
          some_df_1 <- op_lis[[1]]
          some_df_2 <- op_lis[[2]]
          
          #e_l <- e
          
          tf_uni <- unique(some_df_1$from)
          tf_uni <- as.character(tf_uni)
          
          tg_uni <- unique(some_df_1$to)
          tg_uni <- as.character(tg_uni)
          tg_uni[length(tg_uni)+1] <- KO_gene
          
          
        }else{
          ga <- data.frame(stra$arcs)
          
          tg_uni <- c()
          x_l <- which(ga$from %in% KO_gene)
          if(length(x_l) != 0){
            tg_uni <- ga$to[x_l]
          }else{
            tg_uni <- KO_gene
          }
        }
        
        
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l,stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
    }
    
    
    
    
    
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
    
    
    ############################### Lower bounds modification part
    gv <- gpr_vec[[1]]
    xo <- fva_r0_op[[2]]
    
    xo_ind <- which(xo < 0)
    
    for(i in 1:length(xo_ind)){
      xo[xo_ind[i]] <-  xo[xo_ind[i]]*gv[xo_ind[i]]
    }
    
    
    
    # new FVA table - update
    updated_FVA_table <- fva_r0_op
    updated_FVA_table$new_upper_bounds <- new_FVA_ub
    updated_FVA_table$new_lower_bounds <- xo
    
    return(updated_FVA_table)}
  
  
  ########### To check for stoping criteria 
  
  to_check_sink_rxn_iter <- function(fva_pre, fva_1, fva_2,r_begin, r_end,p){
    
    sink_FVA_round_2 <- fva_2[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    
    sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] 
    
    sink_FVA_round_pre <- fva_pre[[3]][r_begin:r_end] 
    
    round_2_vec <- c()
    
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }
    
    round_1_vec <- c()
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }
    
    
    
    # 
    # for(i in 1:length(sink_FVA_round_2)){
    #   if(sink_FVA_round_2[i]> p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 1
    #   }else if(sink_FVA_round_2[i]< p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 0
    #   }
    # }
    
    
    # Maximum fluxes for the sink reactions from round 2 of FVA
    
    # sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    # 
    # round_1_vec <- c()
    # for(i in 1:length(sink_FVA_round_1)){
    #   if(sink_FVA_round_1[i]>0){
    #     round_1_vec[i] = 1
    #   }else if(sink_FVA_round_1[i]<0){
    #     round_1_vec[i] = 0
    #   }
    # }
    
    x_status <- identical(round_1_vec, round_2_vec)
    
    op <- list(x_status, round_1_vec, round_2_vec)
    names(op) <- c("x_status", "Round_1_FVA_max_bin", "Round_2_FVA_max_bin")
    
    return(op)}
  
  
  #######################################################################################################################
  # ########################  Application for the same 
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  library(writexl)
  library(tidyverse)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    
    setwd(curr_wd)
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_3200/"))
    }
    
    grn_SL <- readRDS("Structure_learning.rds")
    
    saveRDS(grn_SL, file = paste0("grn_sl_",count-1,".RDS"))
    
    grn_PL <- readRDS("Parameter_learning.rds")
    
    saveRDS(grn_PL, file = paste0("grn_pl_",count-1,".RDS"))
    
    gge <- read.csv(paste0("Bin_GE_TM2_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    
    setwd(curr_wd)
    matlabr::run_matlab_script("MAT_1b_Initialization_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
    
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM2/grn_3200/"))
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    count <- 1
    
    setwd(curr_wd)
    FVA_to_check <- readxl::read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    colnames(FVA_to_check) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_to_check <- as.data.frame(FVA_to_check)
    
    write.csv(FVA_to_check, paste0("FVA_incorp_",count-1,".csv"))
    
    FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_incorp, paste0("FBA_incorp_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_to_check
    
    FVA_III <- FVA_to_check
    #FVA_III <- FBA_incorp
    
    FVA_II <- FVA_III
    
    xu_n <- xun
    
    
    for(j in 1:length(xu_n)){
      
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_to_check
      
      count <- 1
      
      
      FVA_III <- FVA_to_check
      #FVA_III <- FBA_incorp
      
      FVA_II <- FVA_III
      
      
      FVA_one <- FVA_II
      
      ifva <- Integration_initial(xu_n[j], gene_subsys, grn_SL, grn_PL,gge, curr_wd,count, FVA_up)
      
      
      FVA_P <- ifva 
      #FVA_P <- FBA_incorp
      
      FVA_prev <- FVA_III 
      
      repeat{
        
        #NEW_MGR <- NEW_MGR_prime
        FVA_XP <- FVA_prev
        
        FVA_q <- FVA_P
        
        FVA_i <- FVA_q
        
        NEW_MGR <- get_the_binary_data_iter(met_gene_reg_data,FVA_i,FVA_XP, percen, grn_SL)
        
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(NEW_MGR, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        write.csv(op_intg_1[[1]], paste0("TF_TG_MET_",count+1,".csv"))
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_incorp_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        mup <- read_xlsx(paste0(curr_wd,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
        write.csv(mup, paste0("Updated_FVA_round_incorp_",count+1,".csv"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY_M_3_v2.m", display = TRUE, verbose = TRUE)
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        FVA_iplus1 <- as.data.frame( FVA_iplus1)
        colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        
        #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_incorp_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_incorp_",count+1,".csv"))
        
        
        iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FVA_iplus1,10,11, percen)
        #iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FBA_iplus1,2713,2778, percen)
        
        
        
        write.csv(iter_op$Round_1_FVA_max_bin,paste0("BV_",count,".csv"))
        write.csv(iter_op$Round_2_FVA_max_bin,paste0("BV_",count+1,".csv"))
        
        
        
        count <- count + 1
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        #NEW_MGR_prime <- iter_op$Round_2_FVA_max_bin
        
        FVA_prev <- FVA_i
        
        FVA_P <- FVA_iplus1
        
      }
      
      bv <- c()
      for(i in 0:(count)){
        bv <- c(bv,paste0("BV_",i,".csv"))
      }
      
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_incorp_",i,".csv"))
      }
      
      cc1 <- c()
      for(i in 0:(count)){
        cc1 <- c(cc1,paste0("grn_sl_",i,".RDS"))
      }
      
      cc2 <- c()
      for(i in 0:(count)){
        cc2 <- c(cc2,paste0("grn_pl_",i,".RDS"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_incorp_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_incorp_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      mm <- c()
      for(i in 1:(count)){
        mm <- c(mm,paste0("TF_TG_MET_",i,".csv"))
      }
      
      my_files <- c("FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",vv,bb,pp,uu,mm,cc1,cc2,bv)
      curr_wd_ <- paste0(curr_wd,"/")
      ko_dir <- paste0(curr_wd_,xu_n[j],"/")
      
      dir.create(paste0(curr_wd_,xu_n[j]))
      
      file.rename(from = paste0(curr_wd_, my_files),to = paste0(ko_dir, my_files))
      
    }
    library(ff)
    dir.create(paste0(curr_wd,"/KO_data_TM2_",ee[e]))
    from <- curr_wd            #Current path of your folder
    to   <- curr_wd            #Path you want to move it.
    
    m <- xu_n
    
    for(i in 1:length(m)){
      path1 <- paste0(from,"/",m[i])
      path2 <- paste0(to,"/KO_data_TM2_",ee[e],"/",m[i])
      file.rename(path1,path2)
      
    }
  }
  
}


##################################################################################
### Single KO CF_S

TM3_Single_KO_CF_S <- function(curr_wd,pe,mi, xun){
  
  
  GRN_CS <- function(KO,sl, pl,g, op1){
    
    if(length(op1) == 0 & KO == "WT"){
      yu <- list(sl,pl)
      names(yu) <- c("CS_SL", "CS_PL")
    }else{
      
      if(length(op1) != 0){
        to_d <- op1[[2]]
        
        bh <- vector(mode = "list", length = nrow(to_d))
        names(bh) <- to_d[[1]]
        
        for(i in 1:length(bh)){
          bh[[i]] <- to_d[[3]][i]
        }
        
        if(KO == "WT"){
          bh <- bh
        }else{
          nx <- names(bh)
          if(KO %in% nx){
            
            koi <- which(nx %in% KO)
            bh[[koi]] <- 0
            
          }else{
            bh[[length(bh)+1]] <- 0
            names(bh)[length(bh)] <- KO
          }
          
        }
        
      }else{
        bh <- list()
        bh[[1]] <- 0
        names(bh) <- KO
      }
      
      
      
      #bo <- paste0("mutilated(sl, evidence = list(",noquote(paste(KO)),"= 0))")
      
      #bo_s <- eval(parse(text = bo))
      bo_s <- mutilated(sl,evidence = bh)
      
      bo_p <- bn.fit(bo_s, g, method = "bayes")
      
      if(length(op1) != 0){
        bess <- c(to_d[[1]], KO)
        
        keep_ind <- which(names(bo_p) %in% bess)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        
        
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        
        cptB = matrix(c(0, 1), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_2 <- as.table(cptB)
        
        
        # Replace contents of A with B except at the chosen index
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        
        for(i in 1:length(keep_ind_order)){
          xox <- which(names(bh) %in% keep_ind_order[i])
          
          if(bh[[xox]]==1){
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptB)", sep = "")  
            
            eval(parse(text = so))
          }else{
            so <- paste("bo_p_modified$",keep_ind_order[i]," <- as.table(cptA)", sep = "")  
            
            eval(parse(text = so))
          }
          
        }
        
        
        # bo_p_modified[-keep_ind] <- pl[-keep_ind]
        # 
        # # Modify the content at the unchanged index
        # #bo_p_modified$A <- as.table(cptA)
        # 
        # so <- paste("bo_p_modified$",KO," <- as.table(cptA)", sep = "")  
        # 
        # eval(parse(text = so))
        # 
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
      }else{
        cptA = matrix(c(1, 0), ncol = 2, dimnames = list(NULL, c("0", "1")))
        new_value_1 <- as.table(cptA)
        
        keep_ind <- which(names(bo_p) %in% KO)
        keep_ind_order <- names(bo_p)[keep_ind]
        
        bo_p_modified <- bo_p
        
        
        
        bo_p_modified[-keep_ind] <- pl[-keep_ind]
        
        so <- paste("bo_p_modified$",keep_ind_order," <- as.table(cptA)", sep = "")  
        
        eval(parse(text = so))
        
        yu <- list(bo_s,bo_p_modified)
        names(yu) <- c("CS_SL", "CS_PL")
        
        
      }
      
    }
    
    
    
    
    
    return(yu)}
  
  get_the_binary_data <- function(df1, df2, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>0){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    
    grnarcs <- data.frame(sl$arcs)
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    x <- setdiff(nt[[3]], ft)
    
    xid <- which(nt[[3]] %in% x)
    
    df1 <- nt[-xid,]
    
    
    return(df1)}
  
  # get_the_binary_data_iter(met_gene_reg_data,FVA_i,FVA_XP, percen, grn_SL)
  # df1 <- met_gene_reg_data
  # df2 <- FVA_i
  # df3 <- FVA_XP
  # p <- percen
  # sl <- grn_SL
  
  get_the_binary_data_iter <- function(df1, df2, df3, p, sl){
    bin_vec <- rep(0,times=nrow(df1))
    for(i in 1:nrow(df1)){
      for(j in 1:nrow(df2)){
        if(df1[[5]][i]==df2[[1]][j] && round(df2[[3]][j])>p*df3[[3]][j]){
          bin_vec[i]=1
        }
      }
    }
    df1$bin_data <- bin_vec
    
    grnarcs <- data.frame(sl$arcs)
    
    ft <- union(grnarcs$from, grnarcs$to)
    nt <- df1
    
    #x <- setdiff(nt[[3]], ft)
    #x <- setdiff(ft,nt[[3]])
    #xid <- which(nt[[3]] %in% x)
    xid <- which(nt[[3]] %in% ft)
    
    df1 <- nt[xid,]
    
    
    return(df1)}
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
  # df_1 <- NEW_MGR
  # bn_obj <- grn_SL
  # GS <- gene_subsys
  # KO <- xu_n[j]
  
  Integration_part_1 <- function(df_1, bn_obj, GS, KO){
    all <- df_1[[6]]
    
    if((all(all == 0))== FALSE){
      # removing the unkown interactions
      #new_met_gene_reg_data #df1
      q_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[4]][i]=="?"
        if("?" %in% df_1[[4]][i]){
          q_idx <- c(q_idx,i)
        }
      }
      
      p_idx <- c()
      for(i in 1:nrow(df_1)){
        #df_1[[6]][i]=="0"
        if("0" %in% df_1[[6]][i]){
          p_idx <- c(p_idx,i)
        }
      }
      
      pq_idx <- union(q_idx, p_idx)
      #df_1 <- new_met_gene_reg_data
      if(length(pq_idx)!=0){
        df_1 <- df_1[-pq_idx,]
      }else{
        df_1 <- df_1
      }
      
      
      # to remove those TFs with both "+" and "-" interactions
      # get the unique TFs
      unique_tf <- unique(df_1[[3]])
      
      
      if(length(unique_tf)==length(df_1$`TF-gene`)){
        
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
        
        # more_1_gene <- c()
        # for(i in 1:length(x_id)){
        #   more_1_gene <- c(more_1_gene,names(UNI_TF_name_ids[x_id[i]]))
        # }
        
        # more_1_gene_int <- list()
        # for(i in 1:length(x_id)){
        #   more_1_gene_int[[i]] <- df_1[[4]][UNI_TF_name_ids[[x_id[i]]]]
        # }
        
        # names(more_1_gene_int) <- more_1_gene
        
        # to get a df with only 1 rep of regulation
        only_1_idx <- c()
        for(i in 1:length(UNI_TF_name_ids)){
          if(length(UNI_TF_name_ids[[i]])==1){
            only_1_idx <- c(only_1_idx,i)
          }
        }
        
        req_names <- names(UNI_TF_name_ids)[only_1_idx]
        
        modified_interaction_df <- df_1[which(df_1[[3]] %in% req_names),3:4]
        
        
      } else if (length(unique_tf)!=length(df_1$`TF-gene`)){
        
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
          
          #df_1[[6]][i]==0
          if("0" %in% df_1[[6]][i]){
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
        
        
      }
      
      
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
      tf_tg_df[,ncol(tf_tg_df)+1] <- 0
      colnames(tf_tg_df)[ncol(tf_tg_df)] <- "tg_bigg_sym"
      
      for(i in 1:length(GS$Gene)){
        for(j in 1:length(tf_tg_df$to)){
          if(GS$gene_name[i]==tf_tg_df$to[j]){
            tf_tg_df$tg_bigg_sym[j]=GS$Gene[i]
          }
        }
      }
      
      
      # subset of this modified tf_tg_df - to ensure genes that are not present in the metabolic model are removed  (2c-1c)
      # This is the dataframe to focus on for subsequent analysis (to calculate probabilities)
      tf_tg_df_met_model <- tf_tg_df[tf_tg_df$tg_bigg_sym !=0,]
      
      
      tf_tg_df_met_model[,ncol(tf_tg_df_met_model)+1] <- 0
      colnames(tf_tg_df_met_model)[ncol(tf_tg_df_met_model)] <- "Regulation"
      
      for(i in 1:nrow(modified_interaction_df)){
        for(j in 1:nrow(tf_tg_df_met_model)){
          if(tf_tg_df_met_model$from[j]==modified_interaction_df$`TF-gene`[i]){
            tf_tg_df_met_model$Regulation[j]=modified_interaction_df$Interaction[i]
          }
        }
      }
      
      
      
      
      biis <- c()
      for(i in 1:nrow(modified_interaction_df)){
        if(modified_interaction_df$Interaction[i] == "+"){
          biis[i] <- 1
        }else{
          biis[i] <- 0
        }
      }
      modified_interaction_df$status <- biis
      
      y <- union(bn_arcs$from, bn_arcs$to)
      
      y_ind <- which(modified_interaction_df$`TF-gene` %in% y)
      
      modified_interaction_df <- modified_interaction_df[y_ind,]
      
      
      if(KO == "WT"){
        final_df <- tf_tg_df_met_model
      }else{
        bn_arcs <- bn_obj[[3]]
        bn_arcs <- as.data.frame(bn_arcs)
        
        bn_f <- which(bn_arcs$from %in% KO)
        
        
        if(length(bn_f) != 0){
          
          bn_f_a <- bn_arcs[bn_f,]
          bn_f_a[,ncol(bn_f_a)+1] <- 0
          colnames(bn_f_a)[ncol(bn_f_a)] <- "tg_bigg_sym"
          
          for(i in 1:length(GS$Gene)){
            for(j in 1:length(bn_f_a$to)){
              if(GS$gene_name[i]==bn_f_a$to[j]){
                bn_f_a$tg_bigg_sym[j]=GS$Gene[i]
              }
            }
          }
          
          bn_f_a$Regulation <- rep("-", nrow(bn_f_a))
          
          x0 <- which(bn_f_a[[3]] %in% 0)
          bn_f_a <- bn_f_a[-x0,]
          
          
          final_df <- rbind(tf_tg_df_met_model, bn_f_a)
          
          final_df <- unique(final_df)
        }else{
          final_df <- tf_tg_df_met_model
        }
        
      }
      
      
      op_list <- list(final_df,modified_interaction_df)
      names(op_list) <- c("TF_TG_MET_df","TF_TG_int")
      
    }else{
      op_list <- c()
    }
    
    
    
    return(op_list)}
  
  Integration_initial <- function(KO_gene, GS,stra,para, g ,wds, ct,fu){
    if(KO_gene == "WT"){
      setwd(wds)
      matlabr::run_matlab_script("MAT_1b_Initialization_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      # 
      
      # KO_gene = "gnd"
      # stra = grn_SL
      # para = grn_PL
      # g <- gge
      # 
      
      
    }else{
      
      x <- GRN_CS(KO_gene,stra,para, g,c())
      xcspl <- x$CS_PL
      
      saveRDS(x$CS_SL, file = paste0("grn_sl_",ct,".RDS"))
      saveRDS(x$CS_PL, file = paste0("grn_pl_",ct,".RDS"))
      
      ga <- data.frame(stra$arcs)
      
      tg_g <- c()
      x_l <- which(ga$from %in% KO_gene)
      if(length(x_l) != 0){
        tg_g <- ga$to[x_l]
        tg_g[length(tg_g)+1] <- KO_gene
      }else{
        tg_g <- KO_gene
      }
      
      
      all_tgs_eval_bde <- c()
      #set.seed(61)
      for(i in 1:length(tg_g)){
        ##   should b before ( and after e
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
        #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
        qtxt_1 <- paste0("cpquery(xcspl, ","(",paste(noquote(tg_g[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
        x <- eval(parse(text = qtxt_1))
        all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
      }
      
      
      gene_cp_gg_df <- data.frame(tg_g,all_tgs_eval_bde)
      colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if(gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      # colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      # 
      bid <- c() 
      
      # for(i in 1:nrow(gene_cp_gg_df)){
      #   for(j in 1:nrow(GS)){
      #     if((gene_cp_gg_df$MET_MODEL_TG_genes[i]==GS$gene_name[j]) == TRUE){
      #       #bid <- c(bid, GS$Gene[j])
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(gene_cp_gg_df$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      
      gene_cp_gg_df[is.na(gene_cp_gg_df)] <- 0
      
      int2_res <- list(gene_cp_gg_df)
      names(int2_res) <- c("CP_Final")
      
      
      writexl::write_xlsx(int2_res$CP_Final, paste0(wds,"/CP_round_i.xlsx"),col_names = TRUE)
      
      cp <- read_xlsx(paste0(wds,"/CP_round_i.xlsx"), col_names = FALSE)
      write.csv(cp, paste0("CP_incorp_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
      
      gpr_eval_round_i <- read_xlsx(paste0(wds,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
      gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
      colnames(gpr_eval_round_i) <- "gpr_eval"
      
      # step - 2d - part 2
      new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, fu)
      
      writexl::write_xlsx(new_upd_fva,paste0(wds,"/Updated_FVA_round_i.xlsx"))
      
      mup <- read_xlsx(paste0(wds,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
      write.csv(mup, paste0("Updated_FVA_round_incorp_",ct,".csv"))
      
      
      # step - 2a - run FVA with 0 objective
      setwd(wds)
      matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
      
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      # 
      # FVA_prev <- FVA_III 
      
    }
    return(FVA_incorp)}
  
  #Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
  
  Integration_part_2 <- function(op_lis_int_1, KO_gene, GS,para, strat){
    
    if(KO_gene == "WT" & length(op_lis_int_1)==0){
      vm <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
      colnames(vm) <- c("MET_MODEL_TG_genes", "Probability", "Bigg_symb")
      int2_res <- vector("list",2)
      int2_res[[1]] <- vm
      names(int2_res) <- c("CP_Final","evidence_list")
    }else if(KO_gene == "WT"){
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        some_df_1 <- op_lis[[1]]
        some_df_2 <- op_lis[[2]]
        
        ga <- data.frame(stra$arcs)
        
        # tf_uni <- unique(some_df_1$from)
        # tf_uni <- as.character(tf_uni)
        
        tg_uni <- unique(some_df_1$to)
        tg_uni <- as.character(tg_uni)
        
        
        tf_uni  <- some_df_1$from
        tf_uni <- as.character(tf_uni)
        
        sta_vec <- c()
        for(i in 1:nrow(some_df_1)){
          if(some_df_1$Regulation[i]== "+"){
            sta_vec[i] <- 1
          }else{
            sta_vec[i] <- 0
          }
        }
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(tf_uni[i])),"==",paste(noquote(sta_vec[i])),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l, stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
      
    }else{
      
      # get_the_evi_list <- function(op_lis){
      #   
      #   op_df_1 <- op_lis[[1]]
      #   op_df_2 <- op_lis[[2]]
      #   
      #   req_var <- unique(as.character(op_df_1[[1]]))
      #   extra_var <- op_df_2[[1]]
      #   
      #   needed_var <- intersect(req_var,extra_var)
      #   
      #   symb_vec <- list()
      #   for(j in 1:length(needed_var)){
      #     symb_vec[j] <- paste("1")
      #     
      #   }
      #   names(symb_vec) <- needed_var
      #   
      #   return(symb_vec)}
      # 
      # 
      # e_l <- get_the_evi_list(op_lis_int_1)
      # 
      # 
      # if(length(KO_gene)==0){
      #   e_l <- e_l
      # } else if(KO_gene %in% names(e_l)){
      #   x <- which(names(e_l) %in% KO_gene)
      #   e_l[[x]] <- "0"
      # }else {
      #   e_l <- append(e_l,"0")
      #   names(e_l)[length(e_l)] <- KO_gene
      # }
      
      e_l <- list()
      e_l <- "0"
      names(e_l) <- KO_gene
      
      CP_cal_func <- function(op_lis, e, stra){
        
        if(length(op_lis) != 0){
          some_df_1 <- op_lis[[1]]
          some_df_2 <- op_lis[[2]]
          
          #e_l <- e
          
          tf_uni <- unique(some_df_1$from)
          tf_uni <- as.character(tf_uni)
          
          tg_uni <- unique(some_df_1$to)
          tg_uni <- as.character(tg_uni)
        }else{
          ga <- data.frame(stra$arcs)
          
          tg_uni <- c()
          x_l <- which(ga$from %in% KO_gene)
          if(length(x_l) != 0){
            tg_uni <- ga$to[x_l]
            tg_uni[length(tg_uni)+1] <- KO_gene
          }else{
            tg_uni <- KO_gene
          }
        }
        
        
        
        
        
        all_tgs_eval_bde <- c()
        #set.seed(61)
        for(i in 1:length(tg_uni)){
          ##   should b before ( and after e
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e, method =", paste("'lw'"),")")
          #qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",","evidence = e)")
          qtxt_1 <- paste0("cpquery(para, ","(",paste(noquote(tg_uni[i])),"==",paste("1"),")",",",paste(noquote(KO_gene)),"==",paste("0"),")")
          x <- eval(parse(text = qtxt_1))
          all_tgs_eval_bde <- c(all_tgs_eval_bde,x)
        }
        
        
        gene_cp_gg_df <- data.frame(tg_uni,all_tgs_eval_bde)
        colnames(gene_cp_gg_df) <- c("MET_MODEL_TG_genes","Probability")
        
        
        return(gene_cp_gg_df)}
      
      
      CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l,stra = strat)
      
      # bid <- c() 
      # 
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if(CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j]){
      #       bid <- c(bid, GS$Gene[j])
      #     }
      #   }
      # }
      # 
      # CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      # colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      bid <- c() 
      
      # for(i in 1:nrow(CP_DF_MET_MODEL)){
      #   for(j in 1:nrow(GS)){
      #     if((CP_DF_MET_MODEL$MET_MODEL_TG_genes[i]==GS$gene_name[j])==TRUE){
      #       bid[i] <- GS$Gene[j]
      #     }
      #   }
      # }
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
    }
    
    
    
    
    
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
    
    
    ############################### Lower bounds modification part
    gv <- gpr_vec[[1]]
    xo <- fva_r0_op[[2]]
    
    xo_ind <- which(xo < 0)
    
    for(i in 1:length(xo_ind)){
      xo[xo_ind[i]] <-  xo[xo_ind[i]]*gv[xo_ind[i]]
    }
    
    
    
    # new FVA table - update
    updated_FVA_table <- fva_r0_op
    updated_FVA_table$new_upper_bounds <- new_FVA_ub
    updated_FVA_table$new_lower_bounds <- xo
    
    return(updated_FVA_table)}
  
  
  ########### To check for stoping criteria 
  
  to_check_sink_rxn_iter <- function(fva_pre, fva_1, fva_2,r_begin, r_end,p){
    
    sink_FVA_round_2 <- fva_2[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    
    sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] 
    
    sink_FVA_round_pre <- fva_pre[[3]][r_begin:r_end] 
    
    round_2_vec <- c()
    
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> p*sink_FVA_round_1[i]){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }
    
    round_1_vec <- c()
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> p*sink_FVA_round_pre[i]){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }
    
    
    
    # 
    # for(i in 1:length(sink_FVA_round_2)){
    #   if(sink_FVA_round_2[i]> p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 1
    #   }else if(sink_FVA_round_2[i]< p*sink_FVA_round_1[i]){
    #     round_2_vec[i] = 0
    #   }
    # }
    
    
    # Maximum fluxes for the sink reactions from round 2 of FVA
    
    # sink_FVA_round_1 <- fva_1[[3]][r_begin:r_end] # change the sink reaction indices -- Corresponds to Toy_Model
    # 
    # round_1_vec <- c()
    # for(i in 1:length(sink_FVA_round_1)){
    #   if(sink_FVA_round_1[i]>0){
    #     round_1_vec[i] = 1
    #   }else if(sink_FVA_round_1[i]<0){
    #     round_1_vec[i] = 0
    #   }
    # }
    
    x_status <- identical(round_1_vec, round_2_vec)
    
    op <- list(x_status, round_1_vec, round_2_vec)
    names(op) <- c("x_status", "Round_1_FVA_max_bin", "Round_2_FVA_max_bin")
    
    return(op)}
  
  
  #######################################################################################################################
  # ########################  Application for the same 
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  library(writexl)
  library(tidyverse)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    
    setwd(curr_wd)
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_3200/"))
    }
    
    
    grn_SL <- readRDS("Structure_learning.rds")
    
    saveRDS(grn_SL, file = paste0("grn_sl_",count-1,".RDS"))
    
    grn_PL <- readRDS("Parameter_learning.rds")
    
    saveRDS(grn_PL, file = paste0("grn_pl_",count-1,".RDS"))
    
    gge <- read.csv(paste0("Bin_GE_TM3_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    
    setwd(curr_wd)
    matlabr::run_matlab_script("MAT_1b_Initialization_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
    
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    
    if(ee[e] == 3.2){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_3.2/"))
    }else if(ee[e] == 320){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_320/"))
    }else if(ee[e] == 3200){
      setwd(paste0(curr_wd,"/GRN/TM3/grn_3200/"))
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    count <- 1
    
    setwd(curr_wd)
    FVA_to_check <- readxl::read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    colnames(FVA_to_check) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_to_check <- as.data.frame(FVA_to_check)
    
    write.csv(FVA_to_check, paste0("FVA_incorp_",count-1,".csv"))
    
    FBA_incorp <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_incorp, paste0("FBA_incorp_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_to_check
    
    FVA_III <- FVA_to_check
    #FVA_III <- FBA_incorp
    
    FVA_II <- FVA_III
    
    xu_n <- xun
    
    
    for(j in 1:length(xu_n)){
      
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_to_check
      
      count <- 1
      
      
      FVA_III <- FVA_to_check
      #FVA_III <- FBA_incorp
      
      FVA_II <- FVA_III
      
      
      FVA_one <- FVA_II
      
      ifva <- Integration_initial(xu_n[j], gene_subsys, grn_SL, grn_PL,gge, curr_wd,count, FVA_up)
      
      
      FVA_P <- ifva 
      #FVA_P <- FBA_incorp
      
      FVA_prev <- FVA_III 
      
      repeat{
        
        #NEW_MGR <- NEW_MGR_prime
        FVA_XP <- FVA_prev
        
        FVA_q <- FVA_P
        
        FVA_i <- FVA_q
        
        NEW_MGR <- get_the_binary_data_iter(met_gene_reg_data,FVA_i,FVA_XP, percen, grn_SL)
        
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(NEW_MGR, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        write.csv(op_intg_1[[1]], paste0("TF_TG_MET_",count+1,".csv"))
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          saveRDS(x$CS_SL, file = paste0("grn_sl_",count+1,".RDS"))
          saveRDS(x$CS_PL, file = paste0("grn_pl_",count+1,".RDS"))
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_incorp_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2d_part_1_Iteration_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        mup <- read_xlsx(paste0(curr_wd,"/Updated_FVA_round_i.xlsx"), col_names = FALSE)
        write.csv(mup, paste0("Updated_FVA_round_incorp_",count+1,".csv"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("MAT_2a_FVA_Iteration_TOY_M_4_v2.m", display = TRUE, verbose = TRUE)
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        FVA_iplus1 <- as.data.frame( FVA_iplus1)
        colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        
        #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_incorp_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_incorp_",count+1,".csv"))
        
        
        iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FVA_iplus1,9,9, percen)
        #iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FBA_iplus1,2713,2778, percen)
        
        
        write.csv(iter_op$Round_1_FVA_max_bin,paste0("BV_",count,".csv"))
        write.csv(iter_op$Round_2_FVA_max_bin,paste0("BV_",count+1,".csv"))
        
        
        
        count <- count + 1
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        #NEW_MGR_prime <- iter_op$Round_2_FVA_max_bin
        
        FVA_prev <- FVA_i
        
        FVA_P <- FVA_iplus1
        
      }
      
      bv <- c()
      for(i in 0:(count)){
        bv <- c(bv,paste0("BV_",i,".csv"))
      }
      
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_incorp_",i,".csv"))
      }
      
      cc1 <- c()
      for(i in 0:(count)){
        cc1 <- c(cc1,paste0("grn_sl_",i,".RDS"))
      }
      
      cc2 <- c()
      for(i in 0:(count)){
        cc2 <- c(cc2,paste0("grn_pl_",i,".RDS"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_incorp_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_incorp_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      mm <- c()
      for(i in 1:(count)){
        mm <- c(mm,paste0("TF_TG_MET_",i,".csv"))
      }
      
      my_files <- c("FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",vv,bb,pp,uu,mm,cc1,cc2,bv)
      curr_wd_ <- paste0(curr_wd,"/")
      ko_dir <- paste0(curr_wd_,xu_n[j],"/")
      
      dir.create(paste0(curr_wd_,xu_n[j]))
      
      file.rename(from = paste0(curr_wd_, my_files),to = paste0(ko_dir, my_files))
      
    }
    library(ff)
    dir.create(paste0(curr_wd,"/KO_data_TM3_",ee[e]))
    from <- curr_wd            #Current path of your folder
    to   <- curr_wd            #Path you want to move it.
    
    m <- xu_n
    
    for(i in 1:length(m)){
      path1 <- paste0(from,"/",m[i])
      path2 <- paste0(to,"/KO_data_TM3_",ee[e],"/",m[i])
      file.rename(path1,path2)
      
    }
  }
  
}


#dir.create(paste0(curr_wd,"/P_0"))

ee <- c(3.2, 320, 3200)
s <- c()
for(i in 1:length(ee)){
  s <- c(s,paste0("KO_data_TM1_",ee[i]))
}

ee <- c(3.2, 320, 3200)
s2 <- c()
for(i in 1:length(ee)){
  s2 <- c(s2,paste0("KO_data_TM2_",ee[i]))
}

ee <- c(3.2, 320, 3200)
s3 <- c()
for(i in 1:length(ee)){
  s3 <- c(s3,paste0("KO_data_TM3_",ee[i]))
}

#CF_S_1

KO_vec_1 <- c("WT","B","E","A","X","Z")
TM1_Single_KO_CF_S(curr_wd,0,10,xun = KO_vec_1)  


KO_vec_2 <- c("WT","X","A","I")
TM2_Single_KO_CF_S(curr_wd,0,10,xun = KO_vec_2) 


KO_vec_3 <- c("WT","X","A")
TM3_Single_KO_CF_S(curr_wd,0,10,xun = KO_vec_3)  


###################################################### Together dataframe


##########  Functions
ee <- c(3.2,320,3200)
for(i in 1:length(ee)){
  
  Get_ode_simu_SS1 <- function(sw,condi){
    Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8")
    setwd(sw)
    ODE_s <- read.csv(paste0("ODE_simu_",condi,".csv"), header = TRUE)
    colnames(ODE_s) <- Reactions
    ODE_s <- ODE_s[,1:8]
    return(ODE_s)}
  
  
  
  library(matrixStats)
  library(Metrics)
  
  
  # avg_std_func_ode_simu <- function(df){
  #   
  #   avg_s <- colMeans(df)
  #   std_s <- colSds(as.matrix(df))
  #   
  #   as_df <- data.frame(avg_s,std_s)
  #   colnames(as_df) <- c("Average","SD")
  #   
  #   return(as_df)}
  
  
  Med_mad_func_ode_simu <- function(df){
    
    avg_s <- colMedians(as.matrix(df))
    std_s <- colMads(as.matrix(df))
    
    as_df <- data.frame(avg_s,std_s)
    colnames(as_df) <- c("Median","MAD")
    
    return(as_df)}
  
  
  
  get_avg_or_std_func_SS1 <- function(df1,df2,df3,df4,df5,df6,n){
    
    df1_vec <- c()
    for(i in 1:nrow(df1)){
      df1_vec <- c(df1_vec,df1[i,n])
    }
    df2_vec <- c()
    for(i in 1:nrow(df2)){
      df2_vec <- c(df2_vec,df2[i,n])
    }
    df3_vec <- c()
    for(i in 1:nrow(df3)){
      df3_vec <- c(df3_vec,df3[i,n])
    }
    df4_vec <- c()
    for(i in 1:nrow(df4)){
      df4_vec <- c(df4_vec,df4[i,n])
    }
    df5_vec <- c()
    for(i in 1:nrow(df5)){
      df5_vec <- c(df5_vec,df5[i,n])
    }
    df6_vec <- c()
    for(i in 1:nrow(df6)){
      df6_vec <- c(df6_vec,df6[i,n])
    }
    
    df_vec_tog <- c(df1_vec,df2_vec,df3_vec,df4_vec,df5_vec,df6_vec)
  } 
  
  
  
  pred_output_func_wo_rep_SS1 <- function(wt_df,b_df,e_df,z_df,a_df,x_df,meth){
    
    p_wt <- wt_df[,meth]
    p_b <- b_df[,meth]
    p_e <- e_df[,meth]
    p_z <- z_df[,meth]
    p_a <- a_df[,meth]
    p_x <- x_df[,meth]
    
    all_p <- c(p_wt,p_b,p_e,p_z,p_a,p_z)
    return(all_p)}
  
  
  
  ###############################################
  ############################################### Obtaining the data 
  library(readxl)
  
  
  exch_rate <- ee[i]
  TM = "TM_1"
  
  
  ####################################  CF-S
  # WT
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/WT/"))
  
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_WT_S <- y$V1
  
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/Full_method/based_on_new_GE_data/",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_WT_S <- z[,2]
  FVA_max_WT_S <- z[,3]
  
  # B KO
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/B/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_B_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/B/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_B_S <- z[,2]
  FVA_max_B_S <- z[,3]
  
  
  # E KO
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/E/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_E_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/E/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_E_S <- z[,2]
  FVA_max_E_S <- z[,3]
  
  
  # Z KO
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/Z/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_Z_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/Z/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_Z_S <- z[,2]
  FVA_max_Z_S <- z[,3]
  
  # X KO
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/X/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_X_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/X/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_X_S <- z[,2]
  FVA_max_X_S <- z[,3]
  
  # A KO
  setwd(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/A/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_A_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM1_",exch_rate,"/A/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_A_S <- z[,2]
  FVA_max_A_S <- z[,3]
  

  
  
  ################################################ TRIMER results 
  
  #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Toy_M_1/TRIMER/based_on_new_GE_data/")
  
  setwd(paste0(curr_wd,"Other_dats/TRIMER/TM1/"))
  #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TRIMER/TM1/")
  
  v <- read.csv(paste0("TRIMER_TM_1_",exch_rate,".csv"), header = F)
  
  v <- t(v)
  v <- as.data.frame(v)
  colnames(v) = c("WT","B","E","Z","A","X")
  
  ########################################################  GIMME Results
  
  ####### WT
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/WT/"))
  # setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/WT/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_WT <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_WT <- z[,2]
  GIMME_FVA_max_WT <- z[,3]
  
  ####### A_KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/A_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/A_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_A <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_A <- z[,2]
  GIMME_FVA_max_A <- z[,3]
  
  ####### B_KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/B_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/B_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_B <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_B <- z[,2]
  GIMME_FVA_max_B <- z[,3]
  
  ####### E_KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/E_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/E_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_E <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_E <- z[,2]
  GIMME_FVA_max_E <- z[,3]
  
  ####### Z_KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/Z_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/Z_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_Z <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_Z <- z[,2]
  GIMME_FVA_max_Z <- z[,3]
  
  ####### X_KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM1/",exch_rate,"/X_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_1/GIMME/",exch_rate,"/X_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_X <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_X <- z[,2]
  GIMME_FVA_max_X <- z[,3]
  ###################################################################################
  
  WT_meths_df <- data.frame(FBA_WT_S,FVA_min_WT_S,FVA_max_WT_S,v$WT,GIMME_FBA_WT,GIMME_FVA_min_WT,GIMME_FVA_max_WT)
  colnames(WT_meths_df) <- c("CF_S_FBA_WT","CF_S_FVA_min_WT","CF_S_FVA_max_WT","TRIMER_WT","GIMME_FBA_WT","GIMME_FVA_min_WT","GIMME_FVA_max_WT")
  B_meths_df <- data.frame(FBA_B_S,FVA_min_B_S,FVA_max_B_S,v$B,GIMME_FBA_B,GIMME_FVA_min_B,GIMME_FVA_max_B)
  colnames(B_meths_df) <- c("CF_S_FBA_B","CF_S_FVA_min_B","CF_S_FVA_max_B","TRIMER_B","GIMME_FBA_B","GIMME_FVA_min_B","GIMME_FVA_max_B")
  E_meths_df <- data.frame(FBA_E_S,FVA_min_E_S,FVA_max_E_S,v$E,GIMME_FBA_E,GIMME_FVA_min_E,GIMME_FVA_max_E)
  colnames(E_meths_df) <- c("CF_S_FBA_E","CF_S_FVA_min_E","CF_S_FVA_max_E","TRIMER_E","GIMME_FBA_E","GIMME_FVA_min_E","GIMME_FVA_max_E")
  Z_meths_df <- data.frame(FBA_Z_S,FVA_min_Z_S,FVA_max_Z_S,v$Z,GIMME_FBA_Z,GIMME_FVA_min_Z,GIMME_FVA_max_Z)
  colnames(Z_meths_df) <- c("CF_S_FBA_Z","CF_S_FVA_min_Z","CF_S_FVA_max_Z","TRIMER_Z","GIMME_FBA_Z","GIMME_FVA_min_Z","GIMME_FVA_max_Z")
  A_meths_df <- data.frame(FBA_A_S,FVA_min_A_S,FVA_max_A_S,v$A,GIMME_FBA_A,GIMME_FVA_min_A,GIMME_FVA_max_A)
  colnames(A_meths_df) <- c("CF_S_FBA_A","CF_S_FVA_min_A","CF_S_FVA_max_A","TRIMER_A","GIMME_FBA_A","GIMME_FVA_min_A","GIMME_FVA_max_A")
  X_meths_df <- data.frame(FBA_X_S,FVA_min_X_S,FVA_max_X_S,v$X,GIMME_FBA_X,GIMME_FVA_min_X,GIMME_FVA_max_X)
  colnames(X_meths_df) <- c("CF_S_FBA_X","CF_S_FVA_min_X","CF_S_FVA_max_X","TRIMER_X","GIMME_FBA_X","GIMME_FVA_min_X","GIMME_FVA_max_X")
  
  
  ########## Applying the functions
  
  
  swd <- paste0(curr_wd,"ODE_simu/TM1/",exch_rate,"/")
  library(dplyr)
  
  ODE_simu_WT <- Get_ode_simu_SS1(swd,"WT")
  ODE_simu_WT <- ODE_simu_WT %>% mutate_at(1:8, as.numeric)
  
  
  ODE_simu_B <- Get_ode_simu_SS1(swd,"B")
  ODE_simu_B <- ODE_simu_B %>% mutate_at(1:8, as.numeric)
  
  
  ODE_simu_E <- Get_ode_simu_SS1(swd,"E")
  ODE_simu_E <- ODE_simu_E %>% mutate_at(1:8, as.numeric)
  
  
  ODE_simu_Z <- Get_ode_simu_SS1(swd,"Z")
  ODE_simu_Z <- ODE_simu_Z %>% mutate_at(1:8, as.numeric)
  
  ODE_simu_A <- Get_ode_simu_SS1(swd,"A")
  ODE_simu_A <- ODE_simu_A %>% mutate_at(1:8, as.numeric)
  
  ODE_simu_X <- Get_ode_simu_SS1(swd,"X")
  ODE_simu_X <- ODE_simu_X %>% mutate_at(1:8, as.numeric)
  
  WT_mm_df <- Med_mad_func_ode_simu(ODE_simu_WT)
  B_mm_df <- Med_mad_func_ode_simu(ODE_simu_B)
  E_mm_df <- Med_mad_func_ode_simu(ODE_simu_E)
  Z_mm_df <- Med_mad_func_ode_simu(ODE_simu_Z)
  A_mm_df <- Med_mad_func_ode_simu(ODE_simu_A)
  X_mm_df <- Med_mad_func_ode_simu(ODE_simu_X)
  
  
  all_med_vec <- get_avg_or_std_func_SS1(WT_mm_df,B_mm_df,E_mm_df,Z_mm_df,A_mm_df,X_mm_df,1)
  all_mad_vec <-  get_avg_or_std_func_SS1(WT_mm_df,B_mm_df,E_mm_df,Z_mm_df,A_mm_df,X_mm_df,2)
  
  
  
  pred_Full_fba_S <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,1)
  pred_Full_min_S <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,2)
  pred_Full_max_S <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,3)
  pred_TRIMER <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,4)
  pred_GIMME_fba <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,5)
  pred_GIMME_min <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,6)
  pred_GIMME_max <- pred_output_func_wo_rep_SS1(WT_meths_df,B_meths_df,E_meths_df,Z_meths_df,A_meths_df,X_meths_df,7)
  
  
  
  Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8")
  Reaction_bar_ <- c(rep(Reactions,6))
  
  
  Condition_bar_ <- c(rep("WT",8),rep("B KO",8),rep("E KO",8),rep("Z KO",8),rep("A KO",8),rep("X KO",8))
  
  
  label_min <- rep("CF_S - FVA min",length(pred_Full_min_S))
  label_max <- rep("CF_S - FVA max",length(pred_Full_max_S))
  label_fba <- rep("CF_S - FBA",length(pred_Full_fba_S))
  label_T <- rep("TRIMER",length(pred_TRIMER))
  label_GIMME_min <- rep("GIMME - FVA min",length(pred_GIMME_min))
  label_GIMME_max <- rep("GIMME - FVA max",length(pred_GIMME_max))
  label_GIMME_fba <- rep("GIMME - FBA",length(pred_GIMME_fba))
  
  
  act_pred_min_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_min_S,label_min,Condition_bar_ )
  colnames(act_pred_min_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_max_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_max_S,label_max,Condition_bar_)
  colnames(act_pred_max_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_fba_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_fba_S,label_fba,Condition_bar_)
  colnames(act_pred_fba_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_TRIMER_df <- data.frame(all_med_vec,all_mad_vec,pred_TRIMER,label_T,Condition_bar_)
  colnames(act_pred_TRIMER_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_min_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_min,label_GIMME_min,Condition_bar_ )
  colnames(act_pred_GIMME_min_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_max_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_max,label_GIMME_max,Condition_bar_ )
  colnames(act_pred_GIMME_max_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_fba_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_fba,label_GIMME_fba,Condition_bar_ )
  colnames(act_pred_GIMME_fba_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  
  
  three_df <- rbind(act_pred_min_df_S,act_pred_max_df_S,act_pred_fba_df_S, act_pred_TRIMER_df,act_pred_GIMME_fba_df,act_pred_GIMME_min_df,act_pred_GIMME_max_df)
  
  
  
  library(ggplot2)
  library(ggpubr)
  
  
  TOT_WT <- cbind(WT_meths_df, WT_mm_df$Median)
  colnames(TOT_WT)[8] <- "Actual"
  
  TOT_B <- cbind(B_meths_df, B_mm_df$Median)
  colnames(TOT_B)[8] <- "Actual"
  
  TOT_E <- cbind(E_meths_df, E_mm_df$Median)
  colnames(TOT_E)[8] <- "Actual"
  
  TOT_Z <- cbind(Z_meths_df, Z_mm_df$Median)
  colnames(TOT_Z)[8] <- "Actual"
  
  TOT_A <- cbind(A_meths_df, A_mm_df$Median)
  colnames(TOT_A)[8] <- "Actual"
  
  TOT_X <- cbind(X_meths_df, X_mm_df$Median)
  colnames(TOT_X)[8] <- "Actual"
  
  
  ##################  Together correlation
  library(WRS2)
  
  a <- cor.test(three_df[1:48,1],three_df[1:48,3], method = "spearman")  #CF-S_MIN
  b <- cor.test(three_df[49:96,1],three_df[49:96,3], method = "spearman") #CF-S_MAX
  c <- cor.test(three_df[97:144,1],three_df[97:144,3], method = "spearman") #CF-S_FBA
  
  
  g <- cor.test(three_df[145:192,1],three_df[145:192,3], method = "spearman") #TRIMER
  
  p <- cor.test(three_df[193:240,1],three_df[193:240,3], method = "spearman") #gimme_fba
  q <- cor.test(three_df[241:288,1],three_df[241:288,3], method = "spearman") #gimme_min
  r <- cor.test(three_df[289:336,1],three_df[289:336,3], method = "spearman") #gimme_max
  
  
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,g$estimate,p$estimate,q$estimate,r$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,g$p.value,p$p.value,q$p.value,r$p.value)
  
  h <- data.frame(cor_v, p_v)
  h$cond <- rep("Together", nrow(h))
  h$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  
  ########################  Separate correlation
  
  ################## WT 
  ### FVA
  a <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_min_WT, method = "spearman")
  b <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_max_WT, method = "spearman")
  ### FBA
  c <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FBA_WT, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(WT_mm_df$Median,WT_meths_df$TRIMER_WT, method = "spearman")
  
  e <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FBA_WT, method = "spearman")
  f <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_min_WT, method = "spearman")
  g <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_max_WT, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_WT <- data.frame(cor_v, p_v)
  h_WT$cond <- rep("WT", nrow(h_WT))
  h_WT$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## B KO 
  ### FVA
  a <- cor.test(B_mm_df$Median,B_meths_df$CF_S_FVA_min_B, method = "spearman")
  b <- cor.test(B_mm_df$Median,B_meths_df$CF_S_FVA_max_B, method = "spearman")
  ### FBA
  c <- cor.test(B_mm_df$Median,B_meths_df$CF_S_FBA_B, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(B_mm_df$Median,B_meths_df$TRIMER_B, method = "spearman")
  
  e <- cor.test(B_mm_df$Median,B_meths_df$GIMME_FBA_B, method = "spearman")
  f <- cor.test(B_mm_df$Median,B_meths_df$GIMME_FVA_min_B, method = "spearman")
  g <- cor.test(B_mm_df$Median,B_meths_df$GIMME_FVA_max_B, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_B <- data.frame(cor_v, p_v)
  h_B$cond <- rep("B KO", nrow(h_B))
  h_B$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## E KO 
  ### FVA
  a <- cor.test(E_mm_df$Median,E_meths_df$CF_S_FVA_min_E, method = "spearman")
  b <- cor.test(E_mm_df$Median,E_meths_df$CF_S_FVA_max_E, method = "spearman")
  ### FBA
  c <- cor.test(E_mm_df$Median,E_meths_df$CF_S_FBA_E, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(E_mm_df$Median,E_meths_df$TRIMER_E, method = "spearman")
  
  e <- cor.test(E_mm_df$Median,E_meths_df$GIMME_FBA_E, method = "spearman")
  f <- cor.test(E_mm_df$Median,E_meths_df$GIMME_FVA_min_E, method = "spearman")
  g <- cor.test(E_mm_df$Median,E_meths_df$GIMME_FVA_max_E, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  

  h_E <- data.frame(cor_v, p_v)
  h_E$cond <- rep("E KO", nrow(h_E))
  h_E$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  
  ################## Z KO 
  ### FVA
  a <- cor.test(Z_mm_df$Median,Z_meths_df$CF_S_FVA_min_Z, method = "spearman")
  b <- cor.test(Z_mm_df$Median,Z_meths_df$CF_S_FVA_max_Z, method = "spearman")
  ### FBA
  c <- cor.test(Z_mm_df$Median,Z_meths_df$CF_S_FBA_Z, method = "spearman")
  
  ### TRIMER
  d <- cor.test(Z_mm_df$Median,Z_meths_df$TRIMER_Z, method = "spearman")
  
  e <- cor.test(Z_mm_df$Median,Z_meths_df$GIMME_FBA_Z, method = "spearman")
  f <- cor.test(Z_mm_df$Median,Z_meths_df$GIMME_FVA_min_Z, method = "spearman")
  g <- cor.test(Z_mm_df$Median,Z_meths_df$GIMME_FVA_max_Z, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
 
  h_Z <- data.frame(cor_v, p_v)
  h_Z$cond <- rep("Z KO", nrow(h_Z))
  h_Z$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## X KO 
  ### FVA
  a <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_min_X, method = "spearman")
  b <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_max_X, method = "spearman")
  ### FBA
  c <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FBA_X, method = "spearman")

  ### TRIMER
  d <- cor.test(X_mm_df$Median,X_meths_df$TRIMER_X, method = "spearman")
  
  e <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FBA_X, method = "spearman")
  f <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_min_X, method = "spearman")
  g <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_max_X, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  

  h_X <- data.frame(cor_v, p_v)
  h_X$cond <- rep("X KO", nrow(h_X))
  h_X$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## A KO 
  ### FVA
  a <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_min_A, method = "spearman")
  b <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_max_A, method = "spearman")
  ### FBA
  c <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FBA_A, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(A_mm_df$Median,A_meths_df$TRIMER_A, method = "spearman")
  
  e <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FBA_A, method = "spearman")
  f <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_min_A, method = "spearman")
  g <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_max_A, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
 
  h_A <- data.frame(cor_v, p_v)
  h_A$cond <- rep("A KO", nrow(h_A))
  h_A$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  
  H <- rbind(h,h_WT,h_B,h_E,h_Z,h_A,h_X)
  
  H$ex_o <- rep(exch_rate,nrow(H))
  H$ex_m <- rep(exch_rate,nrow(H))
  
  H$TM <- rep(TM, nrow(H))
  
  dir.create(paste0(curr_wd,"TM1"))
  setwd(paste0(curr_wd,"TM1"))
  write.csv(H, paste0("Cor_PV_RMSE_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
  
  
  #############################################################################
  #############################################################################
  
  ##### saving actual vs pred data 
  
  Actual_all <- c(WT_mm_df$Median, B_mm_df$Median, E_mm_df$Median, Z_mm_df$Median,A_mm_df$Median,X_mm_df$Median)
  
  CF_S_FBA_all <- c(WT_meths_df$CF_S_FBA_WT,B_meths_df$CF_S_FBA_B, E_meths_df$CF_S_FBA_E, Z_meths_df$CF_S_FBA_Z,A_meths_df$CF_S_FBA_A,X_meths_df$CF_S_FBA_X)
  CF_S_FVA_min_all <- c(WT_meths_df$CF_S_FVA_min_WT,B_meths_df$CF_S_FVA_min_B, E_meths_df$CF_S_FVA_min_E, Z_meths_df$CF_S_FVA_min_Z,A_meths_df$CF_S_FVA_min_A,X_meths_df$CF_S_FVA_min_X)
  CF_S_FVA_max_all <- c(WT_meths_df$CF_S_FVA_max_WT,B_meths_df$CF_S_FVA_max_B, E_meths_df$CF_S_FVA_max_E, Z_meths_df$CF_S_FVA_max_Z,A_meths_df$CF_S_FVA_max_A,X_meths_df$CF_S_FVA_max_X)
  
  TRIMER_all <- c(WT_meths_df$TRIMER_WT,B_meths_df$TRIMER_B, E_meths_df$TRIMER_E, Z_meths_df$TRIMER_Z,A_meths_df$TRIMER_A,X_meths_df$TRIMER_X)
  
  GIMME_FBA_all <- c(WT_meths_df$GIMME_FBA_WT,B_meths_df$GIMME_FBA_B, E_meths_df$GIMME_FBA_E, Z_meths_df$GIMME_FBA_Z,A_meths_df$GIMME_FBA_A,X_meths_df$GIMME_FBA_X)
  GIMME_FVA_min_all <- c(WT_meths_df$GIMME_FVA_min_WT,B_meths_df$GIMME_FVA_min_B, E_meths_df$GIMME_FVA_min_E, Z_meths_df$GIMME_FVA_min_Z,A_meths_df$GIMME_FVA_min_A,X_meths_df$GIMME_FVA_min_X)
  GIMME_FVA_max_all <- c(WT_meths_df$GIMME_FVA_max_WT,B_meths_df$GIMME_FVA_max_B, E_meths_df$GIMME_FVA_max_E, Z_meths_df$GIMME_FVA_max_Z,A_meths_df$GIMME_FVA_max_A,X_meths_df$GIMME_FVA_max_X)
  
  AP_dataf <- data.frame(Actual_all, CF_S_FBA_all, CF_S_FVA_min_all,CF_S_FVA_max_all, TRIMER_all,GIMME_FBA_all, GIMME_FVA_min_all,GIMME_FVA_max_all)
  
  AP_dataf$Settings <- Condition_bar_
  
  AP_dataf$ex_o <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$ex_m <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$TM <- rep(TM, nrow(AP_dataf))
  

  write.csv(AP_dataf, paste0("Actual_Pred_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
}


##################### for TM2

ee <- c(3.2,320,3200)
for(i in 1:length(ee)){
  
  
  Get_ode_simu_SS1 <- function(sw,condi){
    Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
    setwd(sw)
    ODE_s <- read.csv(paste0("ODE_simu_",condi,".csv"), header = TRUE)
    colnames(ODE_s) <- Reactions
    ODE_s <- ODE_s[,1:11]
    return(ODE_s)}
  
  
  
  
  library(matrixStats)
  library(Metrics)
  
  
  Med_mad_func_ode_simu <- function(df){
    
    avg_s <- colMedians(as.matrix(df))
    std_s <- colMads(as.matrix(df))
    
    as_df <- data.frame(avg_s,std_s)
    colnames(as_df) <- c("Median","MAD")
    
    return(as_df)}
  
  
  
  
  
  get_avg_or_std_func_SS1 <- function(df1,df2,df3,df4,n){
    
    df1_vec <- c()
    for(i in 1:nrow(df1)){
      df1_vec <- c(df1_vec,df1[i,n])
    }
    df2_vec <- c()
    for(i in 1:nrow(df2)){
      df2_vec <- c(df2_vec,df2[i,n])
    }
    df3_vec <- c()
    for(i in 1:nrow(df3)){
      df3_vec <- c(df3_vec,df3[i,n])
    }
    df4_vec <- c()
    for(i in 1:nrow(df4)){
      df4_vec <- c(df4_vec,df4[i,n])
    }
    df_vec_tog <- c(df1_vec,df2_vec,df3_vec,df4_vec)
  } 
  
  
  
  pred_output_func_wo_rep_SS1 <- function(wt_df,i_df,x_df,a_df,meth){
    
    p_wt <- wt_df[,meth]
    p_i <- i_df[,meth]
    p_x <- x_df[,meth]
    p_a <- a_df[,meth]
    
    all_p <- c(p_wt,p_i,p_x,p_a)
    return(all_p)}
  
  
  
  
  NAs_func <- function(df1, df2, df3){
    row_no <- c()
    
    x_1 <- c()
    x_2 <- c()
    x_3 <- c()
    x_1 <- which(!complete.cases(df1))
    x_2 <- which(!complete.cases(df2))
    x_3 <- which(!complete.cases(df3))
    
    row_no <- c(x_1, x_2,x_3)
    
    df_1_mod <- df1[-row_no,]
    df_2_mod <- df2[-row_no,]
    df_3_mod <- df3[-row_no,]
    
    nas_l <- list(df_1_mod,df_2_mod,df_3_mod)
    names(nas_l) <- c("Simu_WT","Simu_E","Simu_I")
    
    return(nas_l)}
  
  
  
  ###############################################
  ############################################### Obtaining the data 
  library(readxl)
  
  
  exch_rate <- ee[i]
  TM = "TM_2"
  
  
  # WT
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/Full_method/based_on_new_GE_data/",exch_rate,"/WT/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TM2/",exch_rate,"/WT/"))
  setwd(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/WT/"))
  
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_WT_S <- y$V1
  
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/Full_method/based_on_new_GE_data/",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_WT_S <- z[,2]
  FVA_max_WT_S <- z[,3]
  
  # I KO
  setwd(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/I/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_I_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/I/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_I_S <- z[,2]
  FVA_max_I_S <- z[,3]
  
  
  # X KO
  setwd(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/X/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_X_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/X/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_X_S <- z[,2]
  FVA_max_X_S <- z[,3]
  
  # A KO
  setwd(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/A/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_A_S <- y$V1
  
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM2_",exch_rate,"/A/FVA_to_check.xlsx"), col_names  = F)
  z <- as.data.frame(z)
  FVA_min_A_S <- z[,2]
  FVA_max_A_S <- z[,3]
  

  
  ############################################################## TRIMER
  
  # setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TRIMER/TM2/")
  setwd(paste0(curr_wd,"Other_dats/TRIMER/TM2/"))
  
  v <- read.csv(paste0("TRIMER_TM2_",exch_rate,".csv"), header = F)
  #v <- read.csv("TRIMER_TM3_3_2.csv", header = F)
  v <- t(v)
  v <- as.data.frame(v)
  colnames(v) = c("WT","I","X","A")
  
  ########################################################### GIMME Results
  
  ####### WT
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM2/",exch_rate,"/WT/"))
  # setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/GIMME/",exch_rate,"/WT/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_WT <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_WT <- z[,2]
  GIMME_FVA_max_WT <- z[,3]
  
  ####### A KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM2/",exch_rate,"/A_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/GIMME/",exch_rate,"/A_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_A <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_A <- z[,2]
  GIMME_FVA_max_A <- z[,3]
  
  ####### X KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM2/",exch_rate,"/X_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/GIMME/",exch_rate,"/X_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_X <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_X <- z[,2]
  GIMME_FVA_max_X <- z[,3]
  
  ####### I KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM2/",exch_rate,"/I_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/GIMME/",exch_rate,"/I_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_I <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_I <- z[,2]
  GIMME_FVA_max_I <- z[,3]
  
  #############################################################################################
  
  WT_meths_df <- data.frame(FBA_WT_S,FVA_min_WT_S,FVA_max_WT_S,v$WT,GIMME_FBA_WT,GIMME_FVA_min_WT,GIMME_FVA_max_WT)
  colnames(WT_meths_df) <- c("CF_S_FBA_WT","CF_S_FVA_min_WT","CF_S_FVA_max_WT","TRIMER_WT","GIMME_FBA_WT","GIMME_FVA_min_WT","GIMME_FVA_max_WT")
  I_meths_df <- data.frame(FBA_I_S,FVA_min_I_S,FVA_max_I_S,v$I,GIMME_FBA_I,GIMME_FVA_min_I,GIMME_FVA_max_I)
  colnames(I_meths_df) <-  c("CF_S_FBA_I","CF_S_FVA_min_I","CF_S_FVA_max_I","TRIMER_I","GIMME_FBA_I","GIMME_FVA_min_I","GIMME_FVA_max_I")
  X_meths_df <- data.frame(FBA_X_S,FVA_min_X_S,FVA_max_X_S,v$X,GIMME_FBA_X,GIMME_FVA_min_X,GIMME_FVA_max_X)
  colnames(X_meths_df) <- c("CF_S_FBA_X","CF_S_FVA_min_X","CF_S_FVA_max_X","TRIMER_X","GIMME_FBA_X","GIMME_FVA_min_X","GIMME_FVA_max_X")
  A_meths_df <- data.frame(FBA_A_S,FVA_min_A_S,FVA_max_A_S,v$A,GIMME_FBA_A,GIMME_FVA_min_A,GIMME_FVA_max_A)
  colnames(A_meths_df) <- c("CF_S_FBA_A","CF_S_FVA_min_A","CF_S_FVA_max_A","TRIMER_A","GIMME_FBA_A","GIMME_FVA_min_A","GIMME_FVA_max_A")
  
  
  
  ########## Applying the functions
  
  
  swd <- paste0(curr_wd,"ODE_simu/TM2/",exch_rate,"/")
  
  #swd <- paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/ODE_simu/",exch_rate,"/")
  ##### 10000 exchange
  library(dplyr)
  
  ODE_simu_WT <- Get_ode_simu_SS1(swd,"WT")
  ODE_simu_WT <- ODE_simu_WT %>% mutate_at(1:11, as.numeric)
  
  
  ODE_simu_X <- Get_ode_simu_SS1(swd,"X")
  ODE_simu_X <- ODE_simu_X %>% mutate_at(1:11, as.numeric)
  
  
  ODE_simu_I <- Get_ode_simu_SS1(swd,"I")
  ODE_simu_I <- ODE_simu_I %>% mutate_at(1:11, as.numeric)
  
  
  ODE_simu_A <- Get_ode_simu_SS1(swd,"A")
  ODE_simu_A <- ODE_simu_A %>% mutate_at(1:11, as.numeric)
  
  WT_mm_df <- Med_mad_func_ode_simu(ODE_simu_WT)
  X_mm_df <- Med_mad_func_ode_simu(ODE_simu_X)
  I_mm_df <- Med_mad_func_ode_simu(ODE_simu_I)
  A_mm_df <- Med_mad_func_ode_simu(ODE_simu_A)
  
  
  
  all_med_vec <- get_avg_or_std_func_SS1(WT_mm_df,I_mm_df,X_mm_df,A_mm_df,1)
  all_mad_vec <-  get_avg_or_std_func_SS1(WT_mm_df,I_mm_df,X_mm_df,A_mm_df,2)
  
  
  
  pred_Full_fba_S <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,1)
  pred_Full_min_S <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,2)
  pred_Full_max_S <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,3)
  pred_TRIMER <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,4)
  pred_GIMME_fba <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,5)
  pred_GIMME_min <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,6)
  pred_GIMME_max <- pred_output_func_wo_rep_SS1(WT_meths_df,I_meths_df,X_meths_df,A_meths_df,7)
  
  
  
  Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
  Reaction_bar_ <- c(rep(Reactions,4))
  
  
  Condition_bar_ <- c(rep("WT",11),rep("I KO",11),rep("X KO",11),rep("A KO",11))
  
  
  label_min <- rep("CF_S - FVA min",length(pred_Full_min_S))
  label_max <- rep("CF_S - FVA max",length(pred_Full_max_S))
  label_fba <- rep("CF_S - FBA",length(pred_Full_fba_S))
  label_T <- rep("TRIMER",length(pred_TRIMER))
  label_GIMME_min <- rep("GIMME - FVA min",length(pred_GIMME_min))
  label_GIMME_max <- rep("GIMME - FVA max",length(pred_GIMME_max))
  label_GIMME_fba <- rep("GIMME - FBA",length(pred_GIMME_fba))
  
  act_pred_min_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_min_S,label_min,Condition_bar_ )
  colnames(act_pred_min_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_max_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_max_S,label_max,Condition_bar_ )
  colnames(act_pred_max_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_fba_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_fba_S,label_fba,Condition_bar_ )
  colnames(act_pred_fba_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_TRIMER_df <- data.frame(all_med_vec,all_mad_vec,pred_TRIMER,label_T,Condition_bar_ )
  colnames(act_pred_TRIMER_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_min_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_min,label_GIMME_min,Condition_bar_ )
  colnames(act_pred_GIMME_min_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_max_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_max,label_GIMME_max,Condition_bar_ )
  colnames(act_pred_GIMME_max_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_fba_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_fba,label_GIMME_fba,Condition_bar_ )
  colnames(act_pred_GIMME_fba_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  
  three_df <- rbind(act_pred_min_df_S,act_pred_max_df_S ,act_pred_fba_df_S,act_pred_TRIMER_df,act_pred_GIMME_fba_df,act_pred_GIMME_min_df,act_pred_GIMME_max_df)
  
  
  library(ggplot2)
  library(ggpubr)
  
  library(WRS2)
  
  TOT_WT <- cbind(WT_meths_df, WT_mm_df$Median)
  colnames(TOT_WT)[8] <- "Actual"
  
  TOT_I <- cbind(I_meths_df, I_mm_df$Median)
  colnames(TOT_I)[8] <- "Actual"
  
  TOT_X <- cbind(X_meths_df, X_mm_df$Median)
  colnames(TOT_X)[8] <- "Actual"
  
  TOT_A <- cbind(A_meths_df, A_mm_df$Median)
  colnames(TOT_A)[8] <- "Actual"
  
  
  ##################  Together correlation
  library(WRS2)
  
  a <- cor.test(three_df[1:44,1],three_df[1:44,3], method = "spearman") #cfs min
  b <- cor.test(three_df[45:88,1],three_df[45:88,3], method = "spearman") # cfs max
  c <- cor.test(three_df[89:132,1],three_df[89:132,3], method = "spearman") # cfa fba
  g <- cor.test(three_df[133:176,1],three_df[133:176,3], method = "spearman") # trimer
  
  p <- cor.test(three_df[177:220,1],three_df[177:220,3], method = "spearman") #gimme_fba
  q <- cor.test(three_df[221:264,1],three_df[221:264,3], method = "spearman") #gimme_min
  r <- cor.test(three_df[265:308,1],three_df[265:308,3], method = "spearman") #gimme_max
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,g$estimate,p$estimate,q$estimate,r$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,g$p.value,p$p.value,q$p.value,r$p.value)
  
  
  
  h <- data.frame(cor_v, p_v)
  h$cond <- rep("Together", nrow(h))
  h$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ########################  Separate correlation
  
  ################## WT 
  ### FVA
  a <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_min_WT, method = "spearman")
  b <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_max_WT, method = "spearman")
  ### FBA
  c <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FBA_WT, method = "spearman")
  
  ### TRIMER
  d <- cor.test(WT_mm_df$Median,WT_meths_df$TRIMER_WT, method = "spearman")
  
  e <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FBA_WT, method = "spearman")
  f <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_min_WT, method = "spearman")
  g <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_max_WT, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_WT <- data.frame(cor_v, p_v)
  h_WT$cond <- rep("WT", nrow(h_WT))
  h_WT$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## I KO 
  ### FVA
  a <- cor.test(I_mm_df$Median,I_meths_df$CF_S_FVA_min_I, method = "spearman")
  b <- cor.test(I_mm_df$Median,I_meths_df$CF_S_FVA_max_I, method = "spearman")
  ### FBA
  c <- cor.test(I_mm_df$Median,I_meths_df$CF_S_FBA_I, method = "spearman")
  
  ### TRIMER
  d <- cor.test(I_mm_df$Median,I_meths_df$TRIMER_I, method = "spearman")
  
  e <- cor.test(I_mm_df$Median,I_meths_df$GIMME_FBA_I, method = "spearman")
  f <- cor.test(I_mm_df$Median,I_meths_df$GIMME_FVA_min_I, method = "spearman")
  g <- cor.test(I_mm_df$Median,I_meths_df$GIMME_FVA_max_I, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_I <- data.frame(cor_v, p_v)
  h_I$cond <- rep("I KO", nrow(h_I))
  h_I$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## X KO 
  ### FVA
  a <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_min_X, method = "spearman")
  b <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_max_X, method = "spearman")
  ### FBA
  c <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FBA_X, method = "spearman")
  
  ### TRIMER
  d <- cor.test(X_mm_df$Median,X_meths_df$TRIMER_X, method = "spearman")
  
  e <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FBA_X, method = "spearman")
  f <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_min_X, method = "spearman")
  g <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_max_X, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_X <- data.frame(cor_v, p_v)
  h_X$cond <- rep("X KO", nrow(h_X))
  h_X$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## A KO 
  ### FVA
  a <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_min_A, method = "spearman")
  b <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_max_A, method = "spearman")
  ### FBA
  c <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FBA_A, method = "spearman")
  
  ### TRIMER
  d <- cor.test(A_mm_df$Median,A_meths_df$TRIMER_A, method = "spearman")
  
  e <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FBA_A, method = "spearman")
  f <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_min_A, method = "spearman")
  g <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_max_A, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_A <- data.frame(cor_v, p_v)
  h_A$cond <- rep("A KO", nrow(h_A))
  h_A$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  H <- rbind(h,h_WT,h_I,h_X,h_A)
  
  H$ex_o <- rep(exch_rate,nrow(H))
  H$ex_m <- rep(exch_rate,nrow(H))
  
  H$TM <- rep(TM, nrow(H))
  

  dir.create(paste0(curr_wd,"TM2"))
  setwd(paste0(curr_wd,"TM2"))
  
  write.csv(H, paste0("Cor_PV_RMSE_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
  
  #############################################################################
  #############################################################################
  
  ##### saving actual vs pred data 
  
  Actual_all <- c(WT_mm_df$Median, I_mm_df$Median, X_mm_df$Median, A_mm_df$Median)
  
  CF_S_FBA_all <- c(WT_meths_df$CF_S_FBA_WT, I_meths_df$CF_S_FBA_I,X_meths_df$CF_S_FBA_X,A_meths_df$CF_S_FBA_A)
  CF_S_FVA_min_all <- c(WT_meths_df$CF_S_FVA_min_WT,I_meths_df$CF_S_FVA_min_I,X_meths_df$CF_S_FVA_min_XA_meths_df$CF_S_FVA_min_A)
  CF_S_FVA_max_all <- c(WT_meths_df$CF_S_FVA_max_WT,I_meths_df$CF_S_FVA_max_I,X_meths_df$CF_S_FVA_max_X,A_meths_df$CF_S_FVA_max_A)
  
  TRIMER_all <- c(WT_meths_df$TRIMER_WT,I_meths_df$TRIMER_I, X_meths_df$TRIMER_X, A_meths_df$TRIMER_A)
  
  GIMME_FBA_all <- c(WT_meths_df$GIMME_FBA_WT,I_meths_df$GIMME_FBA_I, X_meths_df$GIMME_FBA_X, A_meths_df$GIMME_FBA_A)
  GIMME_FVA_min_all <- c(WT_meths_df$GIMME_FVA_min_WT,I_meths_df$GIMME_FVA_min_I, X_meths_df$GIMME_FVA_min_X, A_meths_df$GIMME_FVA_min_A)
  GIMME_FVA_max_all <- c(WT_meths_df$GIMME_FVA_max_WT,I_meths_df$GIMME_FVA_max_I, X_meths_df$GIMME_FVA_max_X, A_meths_df$GIMME_FVA_max_A)
  
  AP_dataf <- data.frame(Actual_all, CF_S_FBA_all, CF_S_FVA_min_all,CF_S_FVA_max_all, TRIMER_all,GIMME_FBA_all, GIMME_FVA_min_all,GIMME_FVA_max_all)
  
  AP_dataf$Settings <- Condition_bar_
  
  AP_dataf$ex_o <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$ex_m <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$TM <- rep(TM, nrow(AP_dataf))
  
  write.csv(AP_dataf, paste0("Actual_Pred_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
  
}



############### for TM 3

##########  Functions
ee <- c(3.2,320,3200)
for(i in 1:length(ee)){
  
  Get_ode_simu_SS1 <- function(sw,condi){
    Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8","R9")
    setwd(sw)
    ODE_s <- read.csv(paste0("ODE_simu_",condi,".csv"), header = TRUE)
    colnames(ODE_s) <- Reactions
    ODE_s <- ODE_s[,1:9]
    return(ODE_s)}
  
  
  
  library(matrixStats)
  library(Metrics)
  
  
  # avg_std_func_ode_simu <- function(df){
  #   
  #   avg_s <- colMeans(df)
  #   std_s <- colSds(as.matrix(df))
  #   
  #   as_df <- data.frame(avg_s,std_s)
  #   colnames(as_df) <- c("Average","SD")
  #   
  #   return(as_df)}
  
  
  Med_mad_func_ode_simu <- function(df){
    
    avg_s <- colMedians(as.matrix(df))
    std_s <- colMads(as.matrix(df))
    
    as_df <- data.frame(avg_s,std_s)
    colnames(as_df) <- c("Median","MAD")
    
    return(as_df)}
  
  
  
  get_avg_or_std_func_SS1 <- function(df1,df2,df3,n){
    
    df1_vec <- c()
    for(i in 1:nrow(df1)){
      df1_vec <- c(df1_vec,df1[i,n])
    }
    df2_vec <- c()
    for(i in 1:nrow(df2)){
      df2_vec <- c(df2_vec,df2[i,n])
    }
    df3_vec <- c()
    for(i in 1:nrow(df3)){
      df3_vec <- c(df3_vec,df3[i,n])
    }
    
    df_vec_tog <- c(df1_vec,df2_vec,df3_vec)
  } 
  
  
  
  pred_output_func_wo_rep_SS1 <- function(wt_df,a_df,i_df,meth){
    
    p_wt <- wt_df[,meth]
    p_a <- a_df[,meth]
    p_i <- i_df[,meth]
    
    
    all_p <- c(p_wt,p_a,p_i)
    return(all_p)}
  
  
  ###############################################
  ############################################### Obtaining the data 
  library(readxl)
  
  
  exch_rate <- ee[i]
  TM = "TM_3"
  
  
  # WT
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/Full_method/based_on_new_GE_data/",exch_rate,"/WT/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TM3/",exch_rate,"/WT/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/WT/"))
  setwd(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/WT/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_WT_S <- y$V1
  
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_3/New/Full_method/based_on_new_GE_data/",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/WT/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_WT_S <- z[,2]
  FVA_max_WT_S <- z[,3]
  
  
  # X KO
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/X/"))
  setwd(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/X/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_X_S <- y$V1
  
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/X/FVA_to_check.xlsx"), col_names  = F)
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/X/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_X_S <- z[,2]
  FVA_max_X_S <- z[,3]
  
  # A KO
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/A/"))
  setwd(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/A/"))
  
  y <- read.csv("FBA_to_check.csv", header = F)
  FBA_A_S <- y$V1
  
  #z <- read_xlsx(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S/CF_S_TM3/KO_data_TM3_",exch_rate,"/A/FVA_to_check.xlsx"), col_names  = F)
  z <- read_xlsx(paste0(curr_wd,"KO_data_TM3_",exch_rate,"/A/FVA_to_check.xlsx"), col_names  = F)
  
  z <- as.data.frame(z)
  FVA_min_A_S <- z[,2]
  FVA_max_A_S <- z[,3]
  
  
  
  ############################################################ TRIMER 
  #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/CF_pipeline_TMs/TRIMER/TM3/")
  setwd(paste0(curr_wd,"Other_dats/TRIMER/TM3/"))
  
  v <- read.csv(paste0("TRIMER_TM4_",exch_rate,".csv"), header = F)
  v <- t(v)
  v <- as.data.frame(v)
  colnames(v) = c("WT","X","A")
  
  ########################################################  GIMME Results
  
  ####### WT
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/GIMME/",exch_rate,"/WT/"))
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM3/",exch_rate,"/WT/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_WT <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_WT <- z[,2]
  GIMME_FVA_max_WT <- z[,3]
  
  ####### A KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM3/",exch_rate,"/A_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/GIMME/",exch_rate,"/A_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_A <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_A <- z[,2]
  GIMME_FVA_max_A <- z[,3]
  
  ####### X KO
  setwd(paste0(curr_wd,"Other_dats/GIMME/TM3/",exch_rate,"/X_KO/"))
  #setwd(paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/GIMME/",exch_rate,"/X_KO/"))
  y <- read.csv("GIMME_FBA.csv", header = F)
  GIMME_FBA_X <- y$V1
  
  z <- read.csv("GIMME_FVA.csv", header = F)
  GIMME_FVA_min_X <- z[,2]
  GIMME_FVA_max_X <- z[,3]
  
  
  ################################################################
  
  
  WT_meths_df <- data.frame(FBA_WT_S,FVA_min_WT_S,FVA_max_WT_S,v$WT,GIMME_FBA_WT,GIMME_FVA_min_WT,GIMME_FVA_max_WT)
  colnames(WT_meths_df) <- c("CF_S_FBA_WT","CF_S_FVA_min_WT","CF_S_FVA_max_WT","TRIMER_WT","GIMME_FBA_WT","GIMME_FVA_min_WT","GIMME_FVA_max_WT")
  
  X_meths_df <- data.frame(FBA_X_S,FVA_min_X_S,FVA_max_X_S,v$X,GIMME_FBA_X,GIMME_FVA_min_X,GIMME_FVA_max_X)
  colnames(X_meths_df) <- c("CF_S_FBA_X","CF_S_FVA_min_X","CF_S_FVA_max_X","TRIMER_X","GIMME_FBA_X","GIMME_FVA_min_X","GIMME_FVA_max_X")
  A_meths_df <- data.frame(FBA_A_S,FVA_min_A_S,FVA_max_A_S,v$A,GIMME_FBA_A,GIMME_FVA_min_A,GIMME_FVA_max_A)
  colnames(A_meths_df) <- c("CF_S_FBA_A","CF_S_FVA_min_A","CF_S_FVA_max_A","TRIMER_A","GIMME_FBA_A","GIMME_FVA_min_A","GIMME_FVA_max_A")
  
  
  
  
  ########## Applying the functions
  
  
  swd <- paste0(curr_wd,"ODE_simu/TM3/",exch_rate,"/")
  
  #swd <- paste0("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/TOY_M_4/V3/new/ODE_simu/",exch_rate,"/")
  ##### 10000 exchange
  library(dplyr)
  
  ODE_simu_WT <- Get_ode_simu_SS1(swd,"WT")
  ODE_simu_WT <- ODE_simu_WT %>% mutate_at(1:9, as.numeric)
  
  
  ODE_simu_X <- Get_ode_simu_SS1(swd,"X")
  ODE_simu_X <- ODE_simu_X %>% mutate_at(1:9, as.numeric)
  
  
  ODE_simu_A <- Get_ode_simu_SS1(swd,"A")
  ODE_simu_A <- ODE_simu_A %>% mutate_at(1:9, as.numeric)
  
  ######################################################################
  #####################  FOR 320 AND 3200 CASES
  library(DMwR2)
  ODE_simu_WT <- knnImputation(ODE_simu_WT)
  ODE_simu_A <- knnImputation(ODE_simu_A)
  ODE_simu_X <- knnImputation(ODE_simu_X)
  
  
  WT_mm_df <- Med_mad_func_ode_simu(ODE_simu_WT)
  X_mm_df <- Med_mad_func_ode_simu(ODE_simu_X)
  A_mm_df <- Med_mad_func_ode_simu(ODE_simu_A)
  
  
  
  all_med_vec <- get_avg_or_std_func_SS1(WT_mm_df,A_mm_df,X_mm_df,1)
  all_mad_vec <-  get_avg_or_std_func_SS1(WT_mm_df,A_mm_df,X_mm_df,2)
  
  
  
  
  pred_Full_fba_S <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,1)
  pred_Full_min_S <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,2)
  pred_Full_max_S <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,3)
  pred_TRIMER <-  pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,4)
  pred_GIMME_fba <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,5)
  pred_GIMME_min <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,6)
  pred_GIMME_max <- pred_output_func_wo_rep_SS1(WT_meths_df,X_meths_df,A_meths_df,7)
  
  
  
  Reactions <- c("R1","R2","R3","R4","R5","R6","R7","R8","R9")
  Reaction_bar_ <- c(rep(Reactions,3))
  
  
  Condition_bar_ <- c(rep("WT",9),rep("A KO",9),rep("X KO",9))
  
  
  label_min <- rep("CF_S - FVA min",length(pred_Full_min_S))
  label_max <- rep("CF_S - FVA max",length(pred_Full_max_S))
  label_fba <- rep("CF_S - FBA",length(pred_Full_fba_S))
  label_T <- rep("TRIMER",length(pred_TRIMER))
  label_GIMME_min <- rep("GIMME - FVA min",length(pred_GIMME_min))
  label_GIMME_max <- rep("GIMME - FVA max",length(pred_GIMME_max))
  label_GIMME_fba <- rep("GIMME - FBA",length(pred_GIMME_fba))
  
  act_pred_min_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_min_S,label_min,Condition_bar_ )
  colnames(act_pred_min_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_max_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_max_S,label_max,Condition_bar_ )
  colnames(act_pred_max_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_fba_df_S <- data.frame(all_med_vec,all_mad_vec,pred_Full_fba_S,label_fba,Condition_bar_ )
  colnames(act_pred_fba_df_S) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  
  act_pred_TRIMER_df <- data.frame(all_med_vec,all_mad_vec,pred_TRIMER,label_T,Condition_bar_ )
  colnames(act_pred_TRIMER_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_min_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_min,label_GIMME_min,Condition_bar_ )
  colnames(act_pred_GIMME_min_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_max_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_max,label_GIMME_max,Condition_bar_ )
  colnames(act_pred_GIMME_max_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  act_pred_GIMME_fba_df <- data.frame(all_med_vec,all_mad_vec,pred_GIMME_fba,label_GIMME_fba,Condition_bar_ )
  colnames(act_pred_GIMME_fba_df) <- c("Actual","SD","Predicted","Labels","Conditions")
  
  
  
  
  three_df <- rbind(act_pred_min_df_S,act_pred_max_df_S ,act_pred_fba_df_S,act_pred_TRIMER_df,act_pred_GIMME_fba_df,act_pred_GIMME_min_df,act_pred_GIMME_max_df)
  
  
  
  library(ggplot2)
  library(ggpubr)
  
  TOT_WT <- cbind(WT_meths_df, WT_mm_df$Median)
  colnames(TOT_WT)[8] <- "Actual"
  
  
  TOT_X <- cbind(X_meths_df, X_mm_df$Median)
  colnames(TOT_X)[8] <- "Actual"
  
  TOT_A <- cbind(A_meths_df, A_mm_df$Median)
  colnames(TOT_A)[8] <- "Actual"
  
  ##################  Together correlation
  library(WRS2)
  a <- cor.test(three_df[1:27,1],three_df[1:27,3], method = "spearman")
  b <- cor.test(three_df[28:54,1],three_df[28:54,3], method = "spearman")
  c <- cor.test(three_df[55:81,1],three_df[55:81,3], method = "spearman")
  d <- cor.test(three_df[82:108,1],three_df[82:108,3], method = "spearman")
  
  e <- cor.test(three_df[109:135,1],three_df[109:135,3], method = "spearman")
  f <- cor.test(three_df[136:162,1],three_df[136:162,3], method = "spearman")
  g <- cor.test(three_df[163:189,1],three_df[163:189,3], method = "spearman")
  
  # 
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  h <- data.frame(cor_v, p_v)
  h$cond <- rep("Together", nrow(h))
  h$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  
  
  ########################  Separate correlation
  
  ################## WT 
  ### FVA
  a <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_min_WT, method = "spearman")
  b <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FVA_max_WT, method = "spearman")
  ### FBA
  c <- cor.test(WT_mm_df$Median,WT_meths_df$CF_S_FBA_WT, method = "spearman")
  
  ### TRIMER
  d <- cor.test(WT_mm_df$Median,WT_meths_df$TRIMER_WT, method = "spearman")
  
  e <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FBA_WT, method = "spearman")
  f <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_min_WT, method = "spearman")
  g <- cor.test(WT_mm_df$Median,WT_meths_df$GIMME_FVA_max_WT, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
 
  
  h_WT <- data.frame(cor_v, p_v)
  h_WT$cond <- rep("WT", nrow(h_WT))
  h_WT$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## X KO 
  ### FVA
  a <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_min_X, method = "spearman")
  b <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FVA_max_X, method = "spearman")
  ### FBA
  c <- cor.test(X_mm_df$Median,X_meths_df$CF_S_FBA_X, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(X_mm_df$Median,X_meths_df$TRIMER_X, method = "spearman")
  
  e <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FBA_X, method = "spearman")
  f <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_min_X, method = "spearman")
  g <- cor.test(X_mm_df$Median,X_meths_df$GIMME_FVA_max_X, method = "spearman")
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  

  h_X <- data.frame(cor_v, p_v)
  h_X$cond <- rep("X KO", nrow(h_X))
  h_X$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  ################## A KO 
  ### FVA
  a <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_min_A, method = "spearman")
  b <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FVA_max_A, method = "spearman")
  ### FBA
  c <- cor.test(A_mm_df$Median,A_meths_df$CF_S_FBA_A, method = "spearman")
  
  
  ### TRIMER
  d <- cor.test(A_mm_df$Median,A_meths_df$TRIMER_A, method = "spearman")
  
  e <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FBA_A, method = "spearman")
  f <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_min_A, method = "spearman")
  g <- cor.test(A_mm_df$Median,A_meths_df$GIMME_FVA_max_A, method = "spearman")
  
  
  cor_v <- c(a$estimate,b$estimate,c$estimate,d$estimate,e$estimate,f$estimate,g$estimate)
  p_v <- c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
  
  
  h_A <- data.frame(cor_v, p_v)
  h_A$cond <- rep("A KO", nrow(h_A))
  h_A$meths <-   c("CF_S - FVA (min)","CF_S - FVA (max)","CF_S - FBA", 'TRIMER',"GIMME - FBA","GIMME - FVA (min)","GIMME - FVA (max)")
  
  H <- rbind(h,h_WT,h_A,h_X)
  
 
  
  H$ex_o <- rep(exch_rate,nrow(H))
  H$ex_m <- rep(exch_rate,nrow(H))
  
  H$TM <- rep(TM, nrow(H))
  
  dir.create(paste0(curr_wd,"TM3"))
  setwd(paste0(curr_wd,"TM3"))
  
  write.csv(H, paste0("Cor_PV_RMSE_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
  
  
  #############################################################################
  #############################################################################
  
  ##### saving actual vs pred data 
  
  Actual_all <- c(WT_mm_df$Median, X_mm_df$Median, A_mm_df$Median)
  
  CF_S_FBA_all <- c(WT_meths_df$CF_S_FBA_WT,X_meths_df$CF_S_FBA_X,A_meths_df$CF_S_FBA_A)
  CF_S_FVA_min_all <- c(WT_meths_df$CF_S_FVA_min_WT,X_meths_df$CF_S_FVA_min_XA_meths_df$CF_S_FVA_min_A)
  CF_S_FVA_max_all <- c(WT_meths_df$CF_S_FVA_max_WT,X_meths_df$CF_S_FVA_max_X,A_meths_df$CF_S_FVA_max_A)
  
  
  TRIMER_all <- c(WT_meths_df$TRIMER_WT, X_meths_df$TRIMER_X, A_meths_df$TRIMER_A)
  
  GIMME_FBA_all <- c(WT_meths_df$GIMME_FBA_WT, X_meths_df$GIMME_FBA_X, A_meths_df$GIMME_FBA_A)
  GIMME_FVA_min_all <- c(WT_meths_df$GIMME_FVA_min_WT, X_meths_df$GIMME_FVA_min_X, A_meths_df$GIMME_FVA_min_A)
  GIMME_FVA_max_all <- c(WT_meths_df$GIMME_FVA_max_WT, X_meths_df$GIMME_FVA_max_X, A_meths_df$GIMME_FVA_max_A)
  
  AP_dataf <- data.frame(Actual_all, CF_S_FBA_all, CF_S_FVA_min_all,CF_S_FVA_max_all, TRIMER_all,GIMME_FBA_all, GIMME_FVA_min_all,GIMME_FVA_max_all)
  
  AP_dataf$Settings <- Condition_bar_
  
  AP_dataf$ex_o <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$ex_m <- rep(exch_rate,nrow(AP_dataf))
  AP_dataf$TM <- rep(TM, nrow(AP_dataf))
  
  
  setwd(paste0(curr_wd,"TM3"))
  
  write.csv(AP_dataf, paste0("Actual_Pred_data_",exch_rate,"_",TM,".csv"), col.names = T, row.names = F)
  
}



#################################################  Plot set -  1
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
    setwd(paste0(curr_wd,"TM1/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_1.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd(paste0(curr_wd,"TM1/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_1.csv", header = TRUE)
    
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd(paste0(curr_wd,"TM1/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM1")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_1.csv", header = TRUE)
    
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
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
  # ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FBA.pdf")
  # ggsave("CF_MTR_FBA.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MIN.pdf")
  # ggsave("CF_MTR_FVA_MIN.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MAX.pdf")
  # ggsave("CF_MTR_FVA_MAX.jpeg")
  
  
  
}


############################################# 3.2
#sd <- "D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/Results/Sep_act_vs_pred_plot/"

sd <- curr_wd



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
    setwd(paste0(curr_wd,"TM2/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd(paste0(curr_wd,"TM2/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd(paste0(curr_wd,"TM2/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM2/")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_2.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
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
  # ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FBA.pdf")
  # ggsave("CF_MTR_FBA.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MIN.pdf")
  # ggsave("CF_MTR_FVA_MIN.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MAX.pdf")
  # ggsave("CF_MTR_FVA_MAX.jpeg")
  # 
  
  
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
    setwd(paste0(curr_wd,"TM3/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 320){
    setwd(paste0(curr_wd,"TM3/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_320_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }else if(ec == 3200){
    setwd(paste0(curr_wd,"TM3/"))
    #setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_Trial/TM3/")
    act_3.2 <- read.csv("Actual_Pred_data_3200_TM_3.csv", header = TRUE)
    ### 3.2
    d1 <- act_3.2
  }
  
  d1 <- d1[d1$Settings == g,]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
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
  # ggplot(d1, aes(Actual_all, CF_M_FBA_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FBA))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FBA.pdf")
  # ggsave("CF_MTR_FBA.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_min_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA min))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MIN.pdf")
  # ggsave("CF_MTR_FVA_MIN.jpeg")
  # ggplot(d1, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+geom_smooth(method = "lm")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Exchange rate set at ", EO," for the ODE simulations and ", EM, " for \nthe Method under ",COND," condition in ",TM , " (Normalised actual vs predicted)"))+theme_bw()
  # ggsave("CF_MTR_FVA_MAX.pdf")
  # ggsave("CF_MTR_FVA_MAX.jpeg")
  
  
  
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


############################################################################################
######## plot set - 2



cor_dat_df <- function(df1, kv, ex){
  
  # df1$CFMTR_FBA <- df2[[1]]
  # df1$CFMTR_FVA_min <- df2[[2]]
  # df1$CFMTR_FVA_max <- df2[[3]]
  
  x1 <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
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
setwd(paste0(curr_wd,"TM1/"))
###############  3.2
act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_1.csv", header = TRUE)

kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")

a <- cor_dat_df(act_3.2, kv_vec, 3.2)
############### 320
setwd(paste0(curr_wd,"TM1/"))
act_320 <- read.csv("Actual_Pred_data_320_TM_1.csv", header = TRUE)

kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")

b <- cor_dat_df(act_320, kv_vec, 320)

############### 3200
setwd(paste0(curr_wd,"TM1/"))
act_3200 <- read.csv("Actual_Pred_data_3200_TM_1.csv", header = TRUE)

kv_vec <- c("WT","B KO","E KO","Z KO","A KO","X KO")
c <- cor_dat_df(act_3200, kv_vec, 3200)


overall_SC_cfmtr_tm1_per0_corr_data <- rbind(a,b,c)
overall_SC_cfmtr_tm1_per0_corr_data$TMS <- rep("TM1",nrow(overall_SC_cfmtr_tm1_per0_corr_data))

#######################################################################################################
######################   TM2 
setwd(paste0(curr_wd,"TM2/"))
###############  3.2
act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_2.csv", header = TRUE)

kv_vec <- c("WT","I KO","X KO","A KO")

p <- cor_dat_df(act_3.2, kv_vec, 3.2)
############### 320
setwd(paste0(curr_wd,"TM2/"))
act_320 <- read.csv("Actual_Pred_data_320_TM_2.csv", header = TRUE)

kv_vec <- c("WT","I KO","X KO","A KO")

q <- cor_dat_df(act_320, kv_vec, 320)

############### 3200
setwd(paste0(curr_wd,"TM2/"))
act_3200 <- read.csv("Actual_Pred_data_3200_TM_2.csv", header = TRUE)

kv_vec <- c("WT","I KO","X KO","A KO")
r <- cor_dat_df(act_3200, kv_vec, 3200)


overall_SC_cfmtr_tm2_per0_corr_data  <- rbind(p,q,r)
overall_SC_cfmtr_tm2_per0_corr_data$TMS <- rep("TM2",nrow(overall_SC_cfmtr_tm2_per0_corr_data))

#######################################################################################################
######################   TM3 
setwd(paste0(curr_wd,"TM3/"))
###############  3.2
act_3.2 <- read.csv("Actual_Pred_data_3.2_TM_3.csv", header = TRUE)

kv_vec <- c("WT","A KO","X KO")

t <- cor_dat_df(act_3.2, kv_vec, 3.2)
############### 320
setwd(paste0(curr_wd,"TM3/"))
act_320 <- read.csv("Actual_Pred_data_320_TM_3.csv", header = TRUE)

kv_vec <- c("WT","A KO","X KO")

u <- cor_dat_df(act_320, kv_vec, 320)

############### 3200
setwd(paste0(curr_wd,"TM3/"))
act_3200 <- read.csv("Actual_Pred_data_3200_TM_3.csv", header = TRUE)

kv_vec <- c("WT","A KO","X KO")
v <- cor_dat_df(act_3200, kv_vec, 3200)


overall_SC_cfmtr_tm3_per0_corr_data <- rbind(t,u,v)
overall_SC_cfmtr_tm3_per0_corr_data$TMS <- rep("TM3",nrow(overall_SC_cfmtr_tm3_per0_corr_data))




##############################################################################
################################## Calling the functions

SC_TMS_df <- rbind(overall_SC_cfmtr_tm1_per0_corr_data,overall_SC_cfmtr_tm2_per0_corr_data,overall_SC_cfmtr_tm3_per0_corr_data)


library(tidyverse)
library(ggExtra)
library(gridExtra)

library(ggplot2)

somed <- curr_wd


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


# ggplot(SC_TMS_df, aes(x = CF_M_FBA_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-MTR (FBA) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
# ggsave("CF_MTR_FBA_vs_CF_MTR_FVA_MAX.pdf")
# ggsave("CF_MTR_FBA_vs_CF_MTR_FVA_MAX.jpeg")
# 
# ggplot(SC_TMS_df, aes(x = CF_M_FVA_min_all, y = CF_M_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for CF-MTR (FVA min) and Actual")+ylab("Spearman correlation (SC) for CF-MTR (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
# ggsave("CF_MTR_FVA_MIN_vs_CF_MTR_FVA_MAX.pdf")
# ggsave("CF_MTR_FVA_MIN_vs_CF_MTR_FVA_MAX.jpeg")

########################## COMPARING CF-S WITH OTHER METHODS

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("TRIMER_all_vs_CF_S_FVA_MAX.pdf")
ggsave("TRIMER_all_vs_CF_S_FVA_MAX.jpeg")

ggplot(SC_TMS_df, aes(x = TRIMER_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for TRIMER and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05, alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("TRIMER_all_vs_CF_S_FVA_MAX_exch_rates.pdf")
ggsave("TRIMER_all_vs_CF_S_FVA_MAX_exch_rates.jpeg")





ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('purple','aquamarine3', 'darkgoldenrod1'))+labs(color = "Toy Models (TM)")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = TMS),width = 0.05, height = 0.05,alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX.pdf")
ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX.jpeg")

ggplot(SC_TMS_df, aes(x = GIMME_FVA_max_all, y = CF_S_FVA_max_all))+geom_abline()+xlab("Spearman correlation (SC) for GIMME (FVA max) and Actual")+ylab("Spearman correlation (SC) for CF-S (FVA max) and Actual")+ggtitle("Comapring the SC (between methods and actual) across \nall the toy models and conditions")+scale_color_manual(values=c('darkslategray2','deepskyblue3', 'blue'))+labs(color = "Exchange rates")+theme_bw()+theme(legend.position = "bottom", legend.direction  = "horizontal")+geom_jitter(aes(color = as.factor(Exchange)),width = 0.05, height = 0.05,alpha= 1, size = 2)+xlim(-0.5,1.0)+ylim(-0.5,1.0)
ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX_exch_rates.pdf")
ggsave("GIMME_FVA_MAX_vs_CF_S_FVA_MAX_exch_rates.jpeg")




#############################################################################################
############## Plot set 3



all_tms_ex <- function(e){
  t <- c(1, 2, 3)
  di <- data.frame()
  
  for(i in 1:length(t)){
    setwd(paste0(curr_wd,"/TM",t[i],"/"))
    
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


x <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
si <- which(colnames(TMS_EXC_ALL) %in% x)


############## Extracting from TMSE_EXC_ALL
ji <- c(4,12,20,28,36,44,59,70,81,92,100,109,118,
        124,131,139,147,155,163,178,189,200,211,219,228,237,
        242,250,258,266,274,282,297,308,319,330,338,347,356)
TMS_EXC_ALL_ <- TMS_EXC_ALL[ji, ]


x <- c("Actual_all","CF_S_FBA_all","CF_S_FVA_min_all","CF_S_FVA_max_all","TRIMER_all","GIMME_FBA_all","GIMME_FVA_min_all", "GIMME_FVA_max_all" )
si <- which(colnames(TMS_EXC_ALL_) %in% x)


for(i in 1:length(si)){
  TMS_EXC_ALL_[[si[i]]] <- scale(TMS_EXC_ALL_[[si[i]]], center = min(TMS_EXC_ALL_[[si[i]]]), scale = max(TMS_EXC_ALL_[[si[i]]]) - min(TMS_EXC_ALL_[[si[i]]]))
}



somed <- curr_wd


dir.create(paste0(somed,"/BM_ACT_VS_PRED/"))


setwd(paste0(somed,"/BM_ACT_VS_PRED/"))

ggplot(TMS_EXC_ALL_, aes(Actual_all, CF_S_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-S (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("CF_S_FVA_MAX.pdf")
ggsave("CF_S_FVA_MAX.jpeg")

# ggplot(TMS_EXC_ALL_, aes(Actual_all, CF_M_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (CF-MTR (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
# ggsave("CF_MTR_FVA_MAX.pdf")
# ggsave("CF_MTR_FVA_MAX.jpeg")

ggplot(TMS_EXC_ALL_, aes(Actual_all, TRIMER_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (TRIMER)")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("TRIMER.pdf")
ggsave("TRIMER.jpeg")

ggplot(TMS_EXC_ALL_, aes(Actual_all, GIMME_FVA_max_all))+geom_point()+xlab("Actual (Median of ODE simulation)")+ylab("Predicted (GIMME (FVA max))")+stat_cor(method = "spearman", size = 5)+ggtitle(paste0("Normalised actual vs predicted for Biomass reactions \nacross KO cases in various TMS under different exchange rate"))+theme_bw()+ylim(0,1.15)+geom_smooth(method = "lm")
ggsave("GIMME_FVA_max.pdf")
ggsave("GIMME_FVA_max.jpeg")



















