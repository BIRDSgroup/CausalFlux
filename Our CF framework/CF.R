
CF_function <- function(curr_wd,xun,vgval,voval,mi, u, exch_rate, ge,gsl,gpl,mgr,gs,fmr){
  
  GRN_CS <- function(KO,sl, pl,g, op1){
    
    if( 0 %in% length(op1) & "WT" %in% KO){
      yu <- list(sl,pl)
      names(yu) <- c("CS_SL", "CS_PL")
    }else{
      
      if(length(op1) != 0 & is_empty(op1[[1]]) == FALSE){
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
      
      
      bo_s <- mutilated(sl,evidence = bh)
      
      bo_p <- bn.fit(bo_s, g, method = "bayes")
      
      if(length(op1) != 0 & is_empty(op1[[1]]) == FALSE){
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
  
  
  
  
  get_the_binary_data_iter <- function(df1, df2, df3, sl){
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
  
  ####################################### Step - 2c  Integration into GRN  #######################################
  
  ######################################################## This section is divided into two functions
  
  ### Function 1 - to output the functional TF regulation
  
  Integration_part_1 <- function(df_1, bn_obj, GS, KO){
    
    # removing the unkown interactions
    #new_met_gene_reg_data #df1
    q_idx <- c()
    for(i in 1:nrow(df_1)){
      if(df_1[[4]][i]=="?"){
        q_idx <- c(q_idx,i)
      }
    }
    
    p_idx <- c()
    for(i in 1:nrow(df_1)){
      if(df_1[[6]][i]=="0"){
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
    
    
    a <- as.data.frame(bn_obj$arcs)
    xxgenes <- df_1$`TF-gene`
    xgtf <- xxgenes  %in%  a$from
    
    if(all(xgtf == FALSE) == TRUE){
      xre <- xxgenes[xgtf]
      
      df_1 <- df_1[df_1$`TF-gene` %in% xre, ]
    }else{
      df_1 <- df_1
    }
    
    
    if(nrow(df_1) == 0){
      op_list <- list("TF_TG_MET_df" = NULL, "TF_TG_int" = NULL)
      #op_list <- vector("list", 2)
      #names(op_list) <- c("TF_TG_MET_df","TF_TG_int")
    }else{
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
    }
    
    
    
    
    return(op_list)}
  
  
  Integration_initial <- function(KO_gene, GS,stra,para, g ,wds, ct,fu){
    if("WT" %in% KO_gene){
      setwd(wds)
      matlabr::run_matlab_script("M_1b.m", display = TRUE, verbose = TRUE)
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check_P1.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_P1_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check_P1.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_P1_",ct,".csv"))
      
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

      bid <- c()

      b <- match(gene_cp_gg_df$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
      colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
      
      gene_cp_gg_df[is.na(gene_cp_gg_df)] <- 0
      
      int2_res <- list(gene_cp_gg_df)
      names(int2_res) <- c("CP_Final")
      
      
      writexl::write_xlsx(int2_res$CP_Final, paste0(wds,"/CP_round_P1_i.xlsx"),col_names = TRUE)
      
      cp <- read_xlsx(paste0(wds,"/CP_round_P1_i.xlsx"), col_names = FALSE)
      write.csv(cp, paste0("CP_incorp_P1_",ct,".csv"))
      
      
      # step - 2d - part 1
      setwd(wds)
      matlabr::run_matlab_script("M_2d.m", display = TRUE, verbose = TRUE)
      
      gpr_eval_round_i <- read_xlsx(paste0(wds,"/GPR_eval_round_P1_i.xlsx"), col_names = FALSE)
      gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
      colnames(gpr_eval_round_i) <- "gpr_eval"
      
      # step - 2d - part 2
      new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, fu)
      
      writexl::write_xlsx(new_upd_fva,paste0(wds,"/Updated_FVA_round_P1_i.xlsx"))
      
      mup <- read_xlsx(paste0(wds,"/Updated_FVA_round_P1_i.xlsx"), col_names = FALSE)
      write.csv(mup, paste0("Updated_FVA_round_incorp_P1_",ct,".csv"))
      
      
      # step - 2a - run FVA with 0 objective
      setwd(wds)
      matlabr::run_matlab_script("M_2a.m", display = TRUE, verbose = TRUE)
      
      
      FVA_incorp <- read_xlsx(paste0(wds,"/FVA_to_check_P1.xlsx"), col_names = FALSE)
      FVA_incorp <- as.data.frame( FVA_incorp)
      colnames(FVA_incorp) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_incorp, paste0("FVA_incorp_P1_",ct,".csv"))
      
      FBA_incorp <- read_csv("FBA_to_check_P1.csv", col_names = FALSE)
      write.csv(FBA_incorp, paste0("FBA_incorp_P1_",ct,".csv"))
      
      # FVA_P <- FVA_incorp
      # #FVA_P <- FBA_incorp
      #
      # FVA_prev <- FVA_III
      
    }
    return(FVA_incorp)}
  
  
  
  Integration_part_2 <- function(op_lis_int_1, KO_gene, GS,para, strat){
    
    if("WT" %in% KO_gene){
      
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
      
      bid <- c()
      
      b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
      bid <- GS$Gene[b]
      
      CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
      colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
      
      CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
      
      int2_res <- list(CP_DF_MET_MODEL, e_l)
      names(int2_res) <- c("CP_Final","evidence_list")
      
      
    }else{
      
      if(is_empty(op_lis_int_1[[1]]) == TRUE){
        xx <- as.data.frame(strat$arcs)
        
        if(KO_gene %in% xx$from){
          
          xxids <- which(xx$from %in% KO_gene)
          tg_uni <- xx$to[xxids]
          
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
          
          bid <- c()
          
          b <- match(gene_cp_gg_df$MET_MODEL_TG_genes,GS$gene_name)
          bid <- GS$Gene[b]
          
          gene_cp_gg_df[,ncol(gene_cp_gg_df)+1] <- bid
          colnames(gene_cp_gg_df)[ncol(gene_cp_gg_df)] <- "Bigg_symb"
          
          gene_cp_gg_df[is.na(gene_cp_gg_df)] <- 0
          
          int2_res <- list(gene_cp_gg_df)
          names(int2_res) <- c("CP_Final")
          
          
          
        }else{
          
          df1 <- data.frame()
          df2 <- data.frame()
          int2_res <- list(df1,df2)
          names(int2_res) <- c("CP_Final","evidence_list")
          
        }
        
        
      }else{
        get_the_evi_list <- function(op_lis){
          
          op_df_1 <- op_lis[[1]]
          op_df_2 <- op_lis[[2]]
          
          req_var <- unique(as.character(op_df_1[[1]]))
          extra_var <- op_df_2[[1]]
          
          needed_var <- intersect(req_var,extra_var)
          
          symb_vec <- list()
          for(j in 1:length(needed_var)){
            symb_vec[j] <- paste("1")
            
          }
          names(symb_vec) <- needed_var
          
          return(symb_vec)}
        
        
        e_l <- get_the_evi_list(op_lis_int_1)
        
        
        if(length(KO_gene)==0){
          e_l <- e_l
        } else if(KO_gene %in% names(e_l)){
          x <- which(names(e_l) %in% KO_gene)
          e_l[[x]] <- "0"
        }else {
          e_l <- append(e_l,"0")
          names(e_l)[length(e_l)] <- KO_gene
        }
        
        e_l <- list()
        e_l <- "0"
        names(e_l) <- KO_gene
        
        CP_cal_func <- function(op_lis, e){
          
          some_df_1 <- op_lis[[1]]
          some_df_2 <- op_lis[[2]]
          
          #e_l <- e
          
          tf_uni <- unique(some_df_1$from)
          tf_uni <- as.character(tf_uni)
          
          tg_uni <- unique(some_df_1$to)
          tg_uni <- as.character(tg_uni)
          tg_uni[length(tg_uni)+1] <- KO_gene
          
          
          
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
        
        
        CP_DF_MET_MODEL <- CP_cal_func(op_lis_int_1, e_l)
        
        bid <- c()
        
        b <- match(CP_DF_MET_MODEL$MET_MODEL_TG_genes,GS$gene_name)
        bid <- GS$Gene[b]
        
        CP_DF_MET_MODEL[,ncol(CP_DF_MET_MODEL)+1] <- bid
        colnames(CP_DF_MET_MODEL)[ncol(CP_DF_MET_MODEL)] <- "Bigg_symb"
        
        CP_DF_MET_MODEL[is.na(CP_DF_MET_MODEL)] <- 0
        
        int2_res <- list(CP_DF_MET_MODEL, e_l)
        names(int2_res) <- c("CP_Final","evidence_list")
      }
      
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
  
  to_check_sink_rxn_iter <- function(fva_pre, fva_1, fva_2,fr){
    
    sink_FVA_round_2 <- fva_2[[3]][fr] # change the sink reaction indices -- Corresponds to Toy_Model
    
    sink_FVA_round_1 <- fva_1[[3]][fr]
    
    sink_FVA_round_pre <- fva_pre[[3]][fr]
    
    round_2_vec <- c()
    
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> 0){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_2)){
        if(round(sink_FVA_round_2[i])> 0){
          round_2_vec[i] = 1
        }else{
          round_2_vec[i] = 0
        }
      }
    }
    
    round_1_vec <- c()
    
    if(p == 0){
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> 0){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }else{
      for(i in 1:length(sink_FVA_round_1)){
        if(round(sink_FVA_round_1[i])> 0){
          round_1_vec[i] = 1
        }else{
          round_1_vec[i] = 0
        }
      }
    }
    
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
  
  #curr_wd <- c("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/CF_S/")
  
  setwd(curr_wd)
  
  
  vg <-  vgval
  vo <- voval
  maxiter = mi
  
  er <- exch_rate
  write.csv(er, file = "Exch_R.csv", row.names = FALSE)
  write.csv(vg, file = "Exch_G.csv", row.names = FALSE)
  write.csv(vo, file = "Exch_O.csv", row.names = FALSE)
  
  gge <- ge
  grn_SL <- gsl
  grn_PL <- gpl
  met_gene_reg_data <- mgr
  gene_subsys <- gs
  FMR_ <- fmr
  

  setwd(curr_wd)
  matlabr::run_matlab_script("M_1b.m", display = TRUE, verbose = TRUE)
  
  FVA_round_0 <- readxl::read_xlsx(paste0(curr_wd,"/FVA_1b_obj_0_P1.xlsx"), col_names = FALSE)
  colnames(FVA_round_0) <- c("Reactions", "Minimum_flux", "Maximum_flux")
  FVA_round_0 <- as.data.frame(FVA_round_0)
  
  count <- 1

  setwd(curr_wd)
  FVA_to_check <- readxl::read_xlsx(paste0(curr_wd,"/FVA_to_check_P1.xlsx"), col_names = FALSE)
  colnames(FVA_to_check) <- c("Reactions", "Minimum_flux", "Maximum_flux")
  FVA_to_check <- as.data.frame(FVA_to_check)
  
  write.csv(FVA_to_check, paste0("FVA_incorp_P1_",count-1,".csv"))
  
  FBA_incorp <- read_csv("FBA_to_check_P1.csv", col_names = FALSE)
  write.csv(FBA_incorp, paste0("FBA_incorp_P1_",count-1,".csv"))
  
  
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
      
      
      NEW_MGR <- get_the_binary_data_iter(met_gene_reg_data,FVA_i,FVA_XP, grn_SL)
      
      #write.csv(NEW_MGR$bin_data, paste0("Bin_",count,".csv"))
      
      #new_met_gene_reg_data$bin_data <- NEW_MGR
      
      # step - 2c - part 1
      op_intg_1 <- Integration_part_1(NEW_MGR, grn_SL, gene_subsys, xu_n[j])
      # step - 2c - part 2
      
      if(xu_n[j] == "WT"){
        if(is_empty(op_intg_1[[1]]) == TRUE){
          
          CP_Final <- data.frame(0,0,0)
          colnames(CP_Final) <- c("MET_MODEL_TG_genes","Probability","Bigg_symb")
          
          op_intg_2 <- list(CP_Final)
          names(op_intg_2) <- "CP_Final" 
          
        }else{
          x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
          xcspl <- x$CS_PL
          op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys, xcspl, grn_SL)
        }
        
      }else{
        x <- GRN_CS(xu_n[j],grn_SL, grn_PL, gge, op_intg_1)
        xcspl <- x$CS_PL
        op_intg_2 <- Integration_part_2(op_intg_1, KO_gene = xu_n[j], gene_subsys,xcspl, grn_SL) 
      }
      
      
      #write the o/p from 2c to xlsx for matlab function 2
      writexl::write_xlsx(op_intg_2$CP_Final, paste0(curr_wd,"/CP_round_P1_i.xlsx"),col_names = TRUE)
      
      cp <- read_xlsx(paste0(curr_wd,"/CP_round_P1_i.xlsx"), col_names = FALSE)
      write.csv(cp, paste0("CP_incorp_P1_",count+1,".csv"))
      
      
      # step - 2d - part 1
      setwd(curr_wd)
      matlabr::run_matlab_script("M_2d.m", display = TRUE, verbose = TRUE)
      
      gpr_eval_round_i <- read_xlsx(paste0(curr_wd,"/GPR_eval_round_P1_i.xlsx"), col_names = FALSE)
      gpr_eval_round_i <- as.data.frame(gpr_eval_round_i)
      colnames(gpr_eval_round_i) <- "gpr_eval"
      
      # step - 2d - part 2
      new_upd_fva <- to_get_upd_FVA(gpr_eval_round_i, FVA_up)
      
      writexl::write_xlsx(new_upd_fva,paste0(curr_wd,"/Updated_FVA_round_P1_i.xlsx"))
      mup <- read_xlsx(paste0(curr_wd,"/Updated_FVA_round_P1_i.xlsx"), col_names = FALSE)
      write.csv(mup, paste0("Updated_FVA_round_incorp_P1_",count+1,".csv"))
      
      # step - 2a - run FVA with 0 objective
      setwd(curr_wd)
      matlabr::run_matlab_script("M_2a.m", display = TRUE, verbose = TRUE)
      
      
      FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check_P1.xlsx"), col_names = FALSE)
      FVA_iplus1 <- as.data.frame( FVA_iplus1)
      colnames(FVA_iplus1) <- c("Reaction names", "Minimum flux", "Maximum flux")
      
      
      #FVA_iplus1 <- read_xlsx(paste0(curr_wd,"/FVA_to_check.xlsx"), col_names = FALSE)
      write.csv(FVA_iplus1, paste0("FVA_incorp_P1_",count+1,".csv"))
      
      FBA_iplus1 <- read_csv("FBA_to_check_P1.csv", col_names = FALSE)
      write.csv(FBA_iplus1, paste0("FBA_incorp_P1_",count+1,".csv"))
      
      
      iter_op <- to_check_sink_rxn_iter(FVA_XP,FVA_i, FVA_iplus1,FMR_)
    
      count <- count + 1
      
      # step - 2e - break
      if(iter_op[[1]]==TRUE|| count == maxiter){
        break
      }
      
      #NEW_MGR_prime <- iter_op$Round_2_FVA_max_bin
      
      FVA_prev <- FVA_i
      
      FVA_P <- FVA_iplus1
    }
    
    vv <- c()
    for(i in 0:(count)){
      vv <- c(vv,paste0("FVA_incorp_P1_",i,".csv"))
    }
    
    bb <- c()
    for(i in 0:(count)){
      bb <- c(bb,paste0("FBA_incorp_P1_",i,".csv"))
    }
    
    pp <- c()
    for(i in 1:(count)){
      pp <- c(pp,paste0("CP_incorp_P1_",i,".csv"))
    }
    
    uu <- c()
    for(i in 1:(count)){
      uu <- c(uu,paste0("Updated_FVA_round_incorp_P1_",i,".csv"))
    }

    
    my_files <- c("FVA_1b_obj_0_P1.xlsx","CP_round_P1_i.xlsx","FBA_to_check_P1.csv","FVA_to_check_P1.xlsx","GPR_eval_round_P1_i.xlsx","Updated_FVA_round_P1_i.xlsx",vv,bb,pp,uu)
    curr_wd_ <- paste0(curr_wd,"/")
    ko_dir <- paste0(curr_wd_,xu_n[j],"/")
    
    dir.create(paste0(curr_wd_,xu_n[j]))
    
    file.rename(from = paste0(curr_wd_, my_files),to = paste0(ko_dir, my_files))
    
  }
  
  library(ff)
  dir.create(paste0(curr_wd,"/Case_P1_",u))
  from <- curr_wd            #Current path of your folder
  to   <- curr_wd            #Path you want to move it.
  
  m <- xu_n
  
  for(i in 1:length(m)){
    path1 <- paste0(from,"/",m[i])
    path2 <- paste0(to,"/Case_P1_",u,"/",m[i])
    file.rename(path1,path2)
  }
   
}

###################################################################################################

curr_wd <- c("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/Parallel_Runs/CF_S/")


setwd("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/CF_S/GRN_REQ_WO_BIGG/")
GSL <- readRDS("bn_TR.RDS")
GPL <- readRDS("bn_params_TR.RDS")

  ## Metabolic module


load("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/CF_S/GRN_REQ_WO_BIGG/train_data_bn.RData")

GE <- train_bin_set_1
for(i in 1:ncol(GE))
{
  GE[,i] <- as.factor(GE[,i])
}
  ## get all the variables loaded
setwd("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/CF_S/MM_REQ/")
# Our IAF1260 data
MGR <- readRDS("Metabolite_Gene_Regulation.RDS")
# TRIMER ecoli
# met_gene_reg_data <- readRDS("Metabolite_Gene_Regulation_TR.RDS")
GS <- readRDS("New_gene_subsys_iML1515.RDS")
colnames(GS) <- c("Gene","gene_name")


setwd("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/Exchange_rxns/")
FMR <- read.csv("Final_FMR_iml1515.csv")

setwd("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/Causal_Surgery/New_Results/LB/")
RD_check <- read.csv("Overall_Ecoli_data_LB_new_results_with_labels.csv", header = T)

RD_check <- RD_check[67:133,]

tu <- nrow(RD_check)


er_ind <- c()
#fmr
for(i in 1:tu){
  CF_function(curr_wd,RD_check$Gene[i],0,8.5,14.5,10,i,er_ind,ge = GE,gsl = GSL, gpl = GPL, mgr = MGR,gs = GS, fmr = FMR$x)
  
  #CF_S_Version_2_cc(curr_wd,RD_check[[1]][i],0,8.5,14.5,10,i,er_ind,ge = GE,gsl = GSL, gpl = GPL, mgr = MGR,gs = GS, fmr = FMR$x)
    
  #CF_S_Version_2_cc(curr_wd,"WT",0,8.5,14.5,10,i,er_ind,ge = GE,gsl = GSL, gpl = GPL, mgr = MGR,gs = GS, fmr = FMR$x)
  
}



s <- c()
for(i in 1:tu){
  s <- c(s,paste0("Case_P1_",i))
}

dir.create(paste0(curr_wd,"/CF_S_R2"))
from <- curr_wd            #Current path of your folder
to   <- curr_wd            #Path you want to move it.

m <- s

for(i in 1:tu){
   path1 <- paste0(from,"/",m[i])
   path2 <- paste0(to,"/CF_S_R2/",m[i])
  file.rename(path1,path2)
}
