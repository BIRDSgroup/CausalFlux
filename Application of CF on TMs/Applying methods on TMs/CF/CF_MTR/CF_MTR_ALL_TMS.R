TM1_Single_KO_CF_MTR <- function(curr_wd,pe,mi, xun){
  
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
  
  
  
  get_the_binary_data_MTR <- function(r1, df1, df2, per){
    
    bn_vec <- c()
    
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    
    bn_vec_2 <- c()
    
    for(i in 1:nrow(r1)){
      for(j in 1:nrow(bn_df)){
        if(r1[i,2] == bn_df[j,1]){
          bn_vec_2 <- c(bn_vec_2,bn_df[j,2])
        }
      }
    }
    
    r1_new <- r1
    r1_new$bin_data <- bn_vec_2
    
    return(r1_new)} 
  
  
  
  
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
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
    if(KO_gene == "WT"){
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_1b_TM1.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
      
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
      write.csv(cp, paste0("CP_iterations_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_2d_TM1.m", display = TRUE, verbose = TRUE)
      
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
      matlabr::run_matlab_script("CF_MTR_MAT_2a_TM1.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
    }
    return(iter_tr_iplus1)}
  
  
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
  
  og_SC <- function(r1, df1, df2, per){
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    return(bn_df)}
  
  
  is_different <- function(A, B, threshold) {
    if (length(A) != length(B)) {
      stop("Vectors A and B must have the same length")
    }
    
    # Calculate proportion of differences
    diff_proportion <- sum(A != B) / length(A)
    
    # Check if the proportion exceeds the threshold
    return(diff_proportion <= threshold)
  }
  
  
  
  #Iter_SC(iter_tr_some_3,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
  
  
  Iter_SC <- function(df1,df2,og, per,t, ccc){
    
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol,df1[,3] ,df2[,3],bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Max_TR","Iter_TR","Bin_vec")
    
    
    write.csv(bn_df, paste0("Binary_vector_iterations_",ccc,".csv" ))
    
    #x_status <- identical(bn_df[,2], og[,2])
    
    x_status <- is_different(og[,2],bn_df[,4],t) 
    
    
    return(x_status)}
  
  #############################################################################################
  #############################################################  Applying the functions
  #############################################################################################
  
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    xu_n <- xun
    
    
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_3200/")
    }
    
    
    grn_SL <- readRDS("Structure_learning.rds")
    grn_PL <- readRDS("Parameter_learning.rds")
    gge <- read.csv(paste0("Bin_GE_TM1_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    ## Metabolic module
    
    setwd(curr_wd)
    matlabr::run_matlab_script("CF_MTR_MAT_1b_TM1.m", display = TRUE, verbose = TRUE)
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    count <- 1
    
    FVA_bef <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    write.csv(FVA_bef, paste0("FVA_iterations_",count-1,".csv"))
    
    FBA_bef <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_bef, paste0("FBA_iterations_",count-1,".csv"))
    
    
    
    #### Iteration step 
    ## get all the variables loaded
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM1/grn_3200/")
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    
    met_gene_reg_data$Met_symbol <- c("m2[c]", "m4[c]", "m5[c]")
    
    
    ################# set-up for stopping criteria
    setwd(curr_wd )
    init_tr <- read.csv("Initial_Max_Turnover.csv", header = F)
    colnames(init_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    init_tr_main <- init_tr
    
    setwd(curr_wd )
    iter_tr <- read.csv("Iteration_Turnover.csv", header = F)
    colnames(iter_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    iter_tr_some <- iter_tr
    
    #og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_some, percen)
    #og_ec <- alt_sc(init_tr_main,iter_tr_some)
    
    write.csv(init_tr, paste0("TR_init_",count-1,".csv"))
    write.csv(iter_tr_some, paste0("TR_iterations_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_bef
    
    
    iter_tr_some_1 <- iter_tr_some
    
    
    ################## for other KO cases
    
    for(j in 1:length(xu_n)){
      
      iter_tr_some_1 <- iter_tr_some
      # 
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_bef
      
      # 
      count <- 1
      
      
      
      iter_tr_some_2 <- iter_tr_some_1
      
      
      iter_tr_iplus1 <- Integration_initial(xu_n[j],gene_subsys, grn_SL, grn_PL,gge, curr_wd, count,FVA_up)
      
      og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_iplus1, percen)
      og_bn_vec_iter <- og_bn_vec
      
      
      ipp <- iter_tr_iplus1
      
      repeat{
        ipp1 <- ipp
        
        og_bn_vec_iter_iter <- og_bn_vec_iter
        
        #iter_tr_some_3 <- iter_tr_some_2
        
        new_met_gene_reg_data <- get_the_binary_data_MTR(met_gene_reg_data,init_tr_main,ipp1,percen)
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(new_met_gene_reg_data, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = "WT", gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_iterations_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2d_TM1.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2a_TM1.m", display = TRUE, verbose = TRUE)
        
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_iterations_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_iterations_",count+1,".csv"))
        
        # FVA_iplus1 <- as.data.frame( FVA_iplus1)
        # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        setwd(curr_wd )
        iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
        colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
        
        write.csv(iter_tr_iplus1, paste0("TR_iterations_",count+1,".csv"))
        
        
        iter_op <- Iter_SC(init_tr_main,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
        
        og_bn_vec_iter <- read.csv(paste0("Binary_vector_iterations_",count,".csv"))
        og_bn_vec_iter <- og_bn_vec_iter[,c(4:5)]
        
        count <- count + 1
        
        #write.csv(count,"Counts_done.csv", row.names = F)
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        
        
        ipp <- iter_tr_iplus1
      }
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_iterations_",i,".csv"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_iterations_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_iterations_",i,".csv"))
      }
      
      ll <- c()
      for(i in 1:count){
        ll <- c(ll,paste0("Binary_vector_iterations_",i,".csv"))
      }
      
      rr <- c()
      for(i in 0:count){
        rr <- c(rr,paste0("TR_iterations_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      
      my_files <- c("TR_init_0.csv","Initial_Max_Turnover.csv","Iteration_Turnover.csv","FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",ll,vv,bb,pp,rr,uu)
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



TM2_Single_KO_CF_MTR <- function(curr_wd,pe,mi, xun){
  
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
  
  
  
  get_the_binary_data_MTR <- function(r1, df1, df2, per){
    
    bn_vec <- c()
    
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    
    bn_vec_2 <- c()
    
    for(i in 1:nrow(r1)){
      for(j in 1:nrow(bn_df)){
        if(r1[i,2] == bn_df[j,1]){
          bn_vec_2 <- c(bn_vec_2,bn_df[j,2])
        }
      }
    }
    
    r1_new <- r1
    r1_new$bin_data <- bn_vec_2
    
    return(r1_new)} 
  
  
  
  
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
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
    if(KO_gene == "WT"){
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_1b_TM2.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
      
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
      write.csv(cp, paste0("CP_iterations_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_2d_TM2.m", display = TRUE, verbose = TRUE)
      
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
      matlabr::run_matlab_script("CF_MTR_MAT_2a_TM2.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
    }
    return(iter_tr_iplus1)}
  
  
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
  
  og_SC <- function(r1, df1, df2, per){
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    return(bn_df)}
  
  
  is_different <- function(A, B, threshold) {
    if (length(A) != length(B)) {
      stop("Vectors A and B must have the same length")
    }
    
    # Calculate proportion of differences
    diff_proportion <- sum(A != B) / length(A)
    
    # Check if the proportion exceeds the threshold
    return(diff_proportion <= threshold)
  }
  
  
  
  #Iter_SC(iter_tr_some_3,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
  
  
  Iter_SC <- function(df1,df2,og, per,t, ccc){
    
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol,df1[,3] ,df2[,3],bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Max_TR","Iter_TR","Bin_vec")
    
    
    write.csv(bn_df, paste0("Binary_vector_iterations_",ccc,".csv" ))
    
    #x_status <- identical(bn_df[,2], og[,2])
    
    x_status <- is_different(og[,2],bn_df[,4],t) 
    
    
    return(x_status)}
  
  #############################################################################################
  #############################################################  Applying the functions
  #############################################################################################
  
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    xu_n <- xun
    
    
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_3200/")
    }
    
    
    grn_SL <- readRDS("Structure_learning.rds")
    grn_PL <- readRDS("Parameter_learning.rds")
    gge <- read.csv(paste0("Bin_GE_TM3_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    ## Metabolic module
    
    setwd(curr_wd)
    matlabr::run_matlab_script("CF_MTR_MAT_1b_TM2.m", display = TRUE, verbose = TRUE)
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    count <- 1
    
    FVA_bef <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    write.csv(FVA_bef, paste0("FVA_iterations_",count-1,".csv"))
    
    FBA_bef <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_bef, paste0("FBA_iterations_",count-1,".csv"))
    
    
    
    #### Iteration step 
    ## get all the variables loaded
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM2/grn_3200/")
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    
    met_gene_reg_data$Met_symbol <- c("m8[c]", "m8[c]", "m6[c]","m6[c]")
    
    
    ################# set-up for stopping criteria
    setwd(curr_wd )
    init_tr <- read.csv("Initial_Max_Turnover.csv", header = F)
    colnames(init_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    init_tr_main <- init_tr
    
    setwd(curr_wd )
    iter_tr <- read.csv("Iteration_Turnover.csv", header = F)
    colnames(iter_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    iter_tr_some <- iter_tr
    
    #og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_some, percen)
    #og_ec <- alt_sc(init_tr_main,iter_tr_some)
    
    write.csv(init_tr, paste0("TR_init_",count-1,".csv"))
    write.csv(iter_tr_some, paste0("TR_iterations_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_bef
    
    
    iter_tr_some_1 <- iter_tr_some
    
    
    ################## for other KO cases
    
    for(j in 1:length(xu_n)){
      
      iter_tr_some_1 <- iter_tr_some
      # 
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_bef
      
      # 
      count <- 1
      
      
      
      iter_tr_some_2 <- iter_tr_some_1
      
      
      iter_tr_iplus1 <- Integration_initial(xu_n[j],gene_subsys, grn_SL, grn_PL,gge, curr_wd, count,FVA_up)
      
      og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_iplus1, percen)
      og_bn_vec_iter <- og_bn_vec
      
      
      ipp <- iter_tr_iplus1
      
      repeat{
        ipp1 <- ipp
        
        og_bn_vec_iter_iter <- og_bn_vec_iter
        
        #iter_tr_some_3 <- iter_tr_some_2
        
        new_met_gene_reg_data <- get_the_binary_data_MTR(met_gene_reg_data,init_tr_main,ipp1,percen)
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(new_met_gene_reg_data, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = "WT", gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_iterations_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2d_TM2.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2a_TM2.m", display = TRUE, verbose = TRUE)
        
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_iterations_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_iterations_",count+1,".csv"))
        
        # FVA_iplus1 <- as.data.frame( FVA_iplus1)
        # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        setwd(curr_wd )
        iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
        colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
        
        write.csv(iter_tr_iplus1, paste0("TR_iterations_",count+1,".csv"))
        
        
        iter_op <- Iter_SC(init_tr_main,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
        
        og_bn_vec_iter <- read.csv(paste0("Binary_vector_iterations_",count,".csv"))
        og_bn_vec_iter <- og_bn_vec_iter[,c(4:5)]
        
        count <- count + 1
        
        #write.csv(count,"Counts_done.csv", row.names = F)
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        
        
        ipp <- iter_tr_iplus1
      }
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_iterations_",i,".csv"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_iterations_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_iterations_",i,".csv"))
      }
      
      ll <- c()
      for(i in 1:count){
        ll <- c(ll,paste0("Binary_vector_iterations_",i,".csv"))
      }
      
      rr <- c()
      for(i in 0:count){
        rr <- c(rr,paste0("TR_iterations_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      
      my_files <- c("TR_init_0.csv","Initial_Max_Turnover.csv","Iteration_Turnover.csv","FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",ll,vv,bb,pp,rr,uu)
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




TM3_Single_KO_CF_MTR <- function(curr_wd,pe,mi, xun){
  
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
  
  
  
  get_the_binary_data_MTR <- function(r1, df1, df2, per){
    
    bn_vec <- c()
    
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    
    bn_vec_2 <- c()
    
    for(i in 1:nrow(r1)){
      for(j in 1:nrow(bn_df)){
        if(r1[i,2] == bn_df[j,1]){
          bn_vec_2 <- c(bn_vec_2,bn_df[j,2])
        }
      }
    }
    
    r1_new <- r1
    r1_new$bin_data <- bn_vec_2
    
    return(r1_new)} 
  
  
  
  
  
  ####################################### Step - 2c  Integration into GRN  ####################################### 
  
  ######################################################## This section is divided into two functions 
  
  ### Function 1 - to output the functional TF regulation
  
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
    if(KO_gene == "WT"){
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_1b_TM3.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
      
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
      write.csv(cp, paste0("CP_iterations_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("CF_MTR_MAT_2d_TM3.m", display = TRUE, verbose = TRUE)
      
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
      matlabr::run_matlab_script("CF_MTR_MAT_2a_TM3.m", display = TRUE, verbose = TRUE)
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_iterations_",ct,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_iterations_",ct,".csv"))
      
      # FVA_iplus1 <- as.data.frame( FVA_iplus1)
      # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      setwd(wds)
      iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
      colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
      
      write.csv(iter_tr_iplus1, paste0("TR_iterations_",ct,".csv"))
      
    }
    return(iter_tr_iplus1)}
  
  
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
  
  og_SC <- function(r1, df1, df2, per){
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] <- 1
        }else{
          bn_vec[i] <- 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol, bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Bin_vec")
    return(bn_df)}
  
  
  is_different <- function(A, B, threshold) {
    if (length(A) != length(B)) {
      stop("Vectors A and B must have the same length")
    }
    
    # Calculate proportion of differences
    diff_proportion <- sum(A != B) / length(A)
    
    # Check if the proportion exceeds the threshold
    return(diff_proportion <= threshold)
  }
  
  
  
  #Iter_SC(iter_tr_some_3,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
  
  
  Iter_SC <- function(df1,df2,og, per,t, ccc){
    
    bn_vec <- c()
    if(per == 0){
      for(i in 1:nrow(df2)){
        if(df2[i,3] > per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }else{
      for(i in 1:nrow(df2)){
        if(df2[i,3] >= per*(df1[i,3])){
          bn_vec[i] = 1
        }else{
          bn_vec[i] = 0
        }
      }
    }
    
    
    
    bn_df <- data.frame(df2$Met_symbol,df1[,3] ,df2[,3],bn_vec)
    colnames(bn_df) <- c("Met_symbol", "Max_TR","Iter_TR","Bin_vec")
    
    
    write.csv(bn_df, paste0("Binary_vector_iterations_",ccc,".csv" ))
    
    #x_status <- identical(bn_df[,2], og[,2])
    
    x_status <- is_different(og[,2],bn_df[,4],t) 
    
    
    return(x_status)}
  
  #############################################################################################
  #############################################################  Applying the functions
  #############################################################################################
  
  library(readxl)
  library(matlabr)
  library(bnlearn)
  library(dplyr)
  
  maxiter = mi
  
  ee <- c(3.2, 320, 3200)
  
  for(e in 1:length(ee)){
    
    
    setwd(curr_wd)
    write.csv(ee[e], file = "Exch_b.csv", row.names = FALSE)
    
    percen <- pe
    count <- 1
    
    xu_n <- xun
    
    
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_3200/")
    }
    
    
    grn_SL <- readRDS("Structure_learning.rds")
    grn_PL <- readRDS("Parameter_learning.rds")
    gge <- read.csv(paste0("Bin_GE_TM4_",ee[e],".csv"), header = TRUE)
    
    for(i in 1:ncol(gge))
    {
      gge[,i] <- as.factor(gge[,i])
    }
    
    ## Metabolic module
    
    setwd(curr_wd)
    matlabr::run_matlab_script("CF_MTR_MAT_1b_TM3.m", display = TRUE, verbose = TRUE)
    
    FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0.xlsx"), col_names = FALSE)
    colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
    FVA_round_0 <- as.data.frame(FVA_round_0)
    
    count <- 1
    
    FVA_bef <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
    write.csv(FVA_bef, paste0("FVA_iterations_",count-1,".csv"))
    
    FBA_bef <- read_csv("FBA_to_check.csv", col_names = FALSE)
    write.csv(FBA_bef, paste0("FBA_iterations_",count-1,".csv"))
    
    
    
    #### Iteration step 
    ## get all the variables loaded
    if(ee[e] == 3.2){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_3.2/")
    }else if(ee[e] == 320){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_320/")
    }else if(ee[e] == 3200){
      setwd("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR/TM3/grn_3200/")
    }
    
    met_gene_reg_data <- readRDS("met_gene_reg_data.RDS")
    gene_subsys <- readRDS("gene_subsys.RDS")
    
    
    met_gene_reg_data$Met_symbol <- c("m8[c]", "m8[c]", "m8[c]")
    
    
    ################# set-up for stopping criteria
    setwd(curr_wd )
    init_tr <- read.csv("Initial_Max_Turnover.csv", header = F)
    colnames(init_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    init_tr_main <- init_tr
    
    setwd(curr_wd )
    iter_tr <- read.csv("Iteration_Turnover.csv", header = F)
    colnames(iter_tr) <- c("Met_ids","Met_symbol","Max_TR")
    
    iter_tr_some <- iter_tr
    
    #og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_some, percen)
    #og_ec <- alt_sc(init_tr_main,iter_tr_some)
    
    write.csv(init_tr, paste0("TR_init_",count-1,".csv"))
    write.csv(iter_tr_some, paste0("TR_iterations_",count-1,".csv"))
    
    
    FVA_up <- FVA_round_0
    #FVA_up <- FVA_bef
    
    
    iter_tr_some_1 <- iter_tr_some
    
    
    ################## for other KO cases
    
    for(j in 1:length(xu_n)){
      
      iter_tr_some_1 <- iter_tr_some
      # 
      FVA_up <- FVA_round_0
      #FVA_up <- FVA_bef
      
      # 
      count <- 1
      
      
      
      iter_tr_some_2 <- iter_tr_some_1
      
      
      iter_tr_iplus1 <- Integration_initial(xu_n[j],gene_subsys, grn_SL, grn_PL,gge, curr_wd, count,FVA_up)
      
      og_bn_vec <- og_SC(met_gene_reg_data,init_tr_main,iter_tr_iplus1, percen)
      og_bn_vec_iter <- og_bn_vec
      
      
      ipp <- iter_tr_iplus1
      
      repeat{
        ipp1 <- ipp
        
        og_bn_vec_iter_iter <- og_bn_vec_iter
        
        #iter_tr_some_3 <- iter_tr_some_2
        
        new_met_gene_reg_data <- get_the_binary_data_MTR(met_gene_reg_data,init_tr_main,ipp1,percen)
        # step - 2c - part 1
        op_intg_1 <- Integration_part_1(new_met_gene_reg_data, grn_SL, gene_subsys, xu_n[j])
        # step - 2c - part 2
        
        
        if(xu_n[j] == "WT"){
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = "WT", gene_subsys, xcspl, grn_SL)
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL)
          
        }
        
        
        #write the o/p from 2c to xlsx for matlab function 2
        writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_i.xlsx"),col_names = TRUE)
        
        cp <- read_xlsx(paste0(curr_wd,"/CP_round_i.xlsx"), col_names = FALSE)
        write.csv(cp, paste0("CP_iterations_",count+1,".csv"))
        
        
        # step - 2d - part 1
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2d_TM3.m", display = TRUE, verbose = TRUE)
        
        gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_i.xlsx"), col_names = FALSE)
        gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
        colnames(gpr_eval_round_i) <- "gpr_eval"
        
        # step - 2d - part 2
        new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
        
        writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_i.xlsx"))
        
        # step - 2a - run FVA with 0 objective
        setwd(curr_wd)
        matlabr::run_matlab_script("CF_MTR_MAT_2a_TM3.m", display = TRUE, verbose = TRUE)
        
        
        
        FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
        write.csv(FVA_iplus1, paste0("FVA_iterations_",count+1,".csv"))
        
        FBA_iplus1 <- read_csv("FBA_to_check.csv", col_names = FALSE)
        write.csv(FBA_iplus1, paste0("FBA_iterations_",count+1,".csv"))
        
        # FVA_iplus1 <- as.data.frame( FVA_iplus1)
        # colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
        
        setwd(curr_wd )
        iter_tr_iplus1 <- read.csv("Iteration_Turnover.csv", header = F)
        colnames(iter_tr_iplus1) <- c("Met_ids","Met_symbol","Max_TR")
        
        write.csv(iter_tr_iplus1, paste0("TR_iterations_",count+1,".csv"))
        
        
        iter_op <- Iter_SC(init_tr_main,iter_tr_iplus1,og_bn_vec_iter_iter,percen, 0, count)
        
        og_bn_vec_iter <- read.csv(paste0("Binary_vector_iterations_",count,".csv"))
        og_bn_vec_iter <- og_bn_vec_iter[,c(4:5)]
        
        count <- count + 1
        
        #write.csv(count,"Counts_done.csv", row.names = F)
        
        # step - 2e - break
        if(iter_op[[1]]==TRUE|| count == maxiter){
          break
        }
        
        
        
        ipp <- iter_tr_iplus1
      }
      
      vv <- c()
      for(i in 0:(count)){
        vv <- c(vv,paste0("FVA_iterations_",i,".csv"))
      }
      
      bb <- c()
      for(i in 0:(count)){
        bb <- c(bb,paste0("FBA_iterations_",i,".csv"))
      }
      
      pp <- c()
      for(i in 1:(count)){
        pp <- c(pp,paste0("CP_iterations_",i,".csv"))
      }
      
      ll <- c()
      for(i in 1:count){
        ll <- c(ll,paste0("Binary_vector_iterations_",i,".csv"))
      }
      
      rr <- c()
      for(i in 0:count){
        rr <- c(rr,paste0("TR_iterations_",i,".csv"))
      }
      
      uu <- c()
      for(i in 1:(count)){
        uu <- c(uu,paste0("Updated_FVA_round_incorp_",i,".csv"))
      }
      
      
      my_files <- c("TR_init_0.csv","Initial_Max_Turnover.csv","Iteration_Turnover.csv","FVA_1b_obj_0.xlsx","CP_round_i.xlsx","FBA_to_check.csv","FVA_to_check.xlsx","GPR_eval_round_i.xlsx","Updated_FVA_round_i.xlsx",ll,vv,bb,pp,rr,uu)
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

library(tidyverse)
curr_wd <- c("D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR_1/")

## Choose a binarization percentage from (0,0.25,0.5,0.75,1)
p <- 0.5 
## Give your value for maximum iterations or keep at 5 used in this work
maxi <- 5

KO_vec_1 <- c("WT","B","E","A","X","Z")
TM1_Single_KO_CF_MTR(curr_wd,p,maxi,xun = KO_vec_1)  


KO_vec_2 <- c("WT","X","A","I")
TM2_Single_KO_CF_MTR(curr_wd,p,maxi,xun = KO_vec_2) 


KO_vec_3 <- c("WT","X","A")
TM3_Single_KO_CF_MTR(curr_wd,p,maxi,xun = KO_vec_3)  

dir.create(paste0(curr_wd,"/Res"))
from <- curr_wd            #Current path of your folder
to   <- curr_wd            #Path you want to move it.

m <- c(s,s2,s3)
tu <- length(m)


for(i in 1:tu){
  path1 <- paste0(from,"/",m[i])
  path2 <- paste0(to,"/Res/",m[i])
  file.rename(path1,path2)
}



















