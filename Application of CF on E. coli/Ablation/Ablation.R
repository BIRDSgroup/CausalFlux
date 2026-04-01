##############################################################################################
#######  pre

## get all the variables loaded
setwd("D:/work/Integrated_network_model/Git_hub_codes/Ecoli/Ablation/")
MGR <- readRDS("Metabolite_Gene_Regulation.RDS")
genes_mgr <- unique(MGR$`TF-gene`)  #57
GSL <- readRDS("bn_TR.RDS")
GRN_arcs <- data.frame(GSL$arcs)
Ecoli_Abla_df <- read.csv("Ecoli_Ablation_Predictions.csv", row.names = 1)

############################################################################################################################ 
############################################################## Shortest path analysis
############################################################################################################################ 

library(igraph)

g <- graph_from_data_frame(GRN_arcs, directed = FALSE)


V2 <- Ecoli_Abla_df$Gene



setdiff(genes_mgr, V(g)$name)   # nodes in V1 not present in graph ## "dsdC" "yqjI" "alaS"
setdiff(V2, V(g)$name)   # nodes in V2 not present in graph  ## 0 

genes_mgr_mod_ind <- which(genes_mgr %in% c("dsdC", "yqjI" ,"alaS"))   ###these genes of MGR not present in GRN
genes_mgr_mod <- genes_mgr[-genes_mgr_mod_ind]

library(igraph)

g <- graph_from_data_frame(GRN_arcs, directed = FALSE)


# compute distances from all V2 nodes (rows) to all V1 nodes (columns)
dist_matrix <- distances(g, v = V2, to = genes_mgr_mod, mode = "out")  # or "all"/"in" depending on direction


# shortest distance from each V2 node to its nearest V1 node
min_dist_to_V1 <- apply(dist_matrix, 1, min, na.rm = TRUE)

# store results in a data frame
result_df <- data.frame(Node = V2, ShortestDistToV1 = min_dist_to_V1)




################## Cols sums of dist_matrix
dist_matrix_1 <- dist_matrix
dist_matrix_1 <- as.data.frame(dist_matrix_1)
dist_matrix_1[dist_matrix_1 == Inf] <- 10



########################################################################################
############ Ablation plots

####################### crp

neigh_crp <- rownames(dist_matrix_1[dist_matrix_1$crp == 1,])
crp_ids <- match(neigh_crp,Ecoli_Abla_df$Gene)
crp_shp <- ifelse(1:nrow(Ecoli_Abla_df) %in% crp_ids, "1_hop", "Normal")
Ecoli_Abla_df$hop_1_crp <- crp_shp



labs_req <- Ecoli_Abla_df[Ecoli_Abla_df$CF_final_labels == 1 & Ecoli_Abla_df$Abla_crp_BM_Final_labs ==0,1] #"astC" "gabT" "cysA" "serC" "cysW" "cysP" "cysU" "argD"
highlight_ids <- match(labs_req, Ecoli_Abla_df$Gene)


intersect(labs_req,neigh_crp)  ## "gabT" "serC"

library(ggrepel)
library(ggpubr)

pos <- position_jitter(width = 0.1, height = 0.1)


ggplot(Ecoli_Abla_df, aes(CF_S_pred, Abla_crp_BM_Final))+geom_jitter(position = pos,aes(color = as.factor(True_labels), shape = as.factor(hop_1_crp)), size = 2, alpha = 0.4)+geom_abline()+theme_bw()+
  labs(title = "Biomass of gene KOs before vs after feedback to crp gene is removed",x = "BM before feedback to crp is removed", y = "BM after feedback to crp is removed")+
  scale_color_manual(values = c("1" = "turquoise4", "0" = "firebrick2"))+
  scale_shape_manual(values=c(3, 16))+
  geom_text_repel(
    data = subset(Ecoli_Abla_df, Gene %in% labs_req),
    aes(label = Gene),
    vjust = -1, # adjust vertical position
    position = pos,
    color = "darkmagenta", # optional highlight color
    min.segment.length = 0,segment.color = NA
  )+stat_cor(method="spearman")

ggsave("Abla_crp.pdf")
ggsave("Abla_crp.jpeg")


##################### fur

neigh_fur <- rownames(dist_matrix_1[dist_matrix_1$fur == 1,])
fur_ids <- match(neigh_fur,Ecoli_Abla_df$Gene)
fur_shp <- ifelse(1:nrow(Ecoli_Abla_df) %in% fur_ids, "1_hop", "Normal")
Ecoli_Abla_df$hop_1_fur <- fur_shp



labs_req <- Ecoli_Abla_df[Ecoli_Abla_df$CF_final_labels ==1 & Ecoli_Abla_df$Abla_fur_BM_Final_labs ==0,1]
#ecoli_check[ecoli_check$CF_S_pred != 0 & ecoli_check$Perturbed_MGR_fur_BM_Final ==0,2]
highlight_ids <- match(labs_req, Ecoli_Abla_df$Gene)


intersect(labs_req,neigh_fur)  ## 0 genes

library(ggrepel)

pos <- position_jitter(width = 0.1, height = 0.1)


ggplot(Ecoli_Abla_df, aes(CF_S_pred, Abla_fur_BM_Final))+geom_jitter(position = pos,aes(color = as.factor(True_labels), shape = as.factor(hop_1_fur)), size = 2, alpha = 0.4)+geom_abline()+theme_bw()+
  labs(title = "Biomass of gene KOs before vs after feedback to fur gene is removed",x = "BM before feedback to fur is removed", y = "BM after feedback to fur is removed")+
  scale_color_manual(values = c("1" = "turquoise4", "0" = "firebrick2"))+
  scale_shape_manual(values=c(3, 16))+
  geom_text_repel(
    data = subset(Ecoli_Abla_df, Gene %in% labs_req),
    aes(label = Gene),
    vjust = -1, # adjust vertical position
    position = pos,
    color = "green", # optional highlight color
    min.segment.length = 0,segment.color = NA
  )+stat_cor(method="spearman")


ggsave("Abla_fur.pdf")
ggsave("Abla_fur.jpeg")



######################## cra

neigh_cra <- rownames(dist_matrix_1[dist_matrix_1$cra == 1,])
cra_ids <- match(neigh_cra,Ecoli_Abla_df$Gene)
cra_shp <- ifelse(1:nrow(Ecoli_Abla_df) %in% cra_ids, "1_hop", "Normal")
Ecoli_Abla_df$hop_1_cra <- cra_shp



labs_req <- Ecoli_Abla_df[Ecoli_Abla_df$CF_final_labels ==1 & Ecoli_Abla_df$Abla_cra_BM_Final_labs ==0,1]
highlight_ids <- match(labs_req, Ecoli_Abla_df$Gene)


intersect(labs_req,neigh_cra)  ## "gabT" "serC"

library(ggrepel)

pos <- position_jitter(width = 0.1, height = 0.1)


ggplot(Ecoli_Abla_df, aes(CF_S_pred, Abla_cra_BM_Final))+geom_jitter(position = pos,aes(color = as.factor(True_labels), shape = as.factor(hop_1_cra)), size = 2, alpha = 0.4)+geom_abline()+theme_bw()+
  labs(title = "Biomass of gene KOs before vs after feedback to cra gene is removed",x = "BM before feedback to cra is removed", y = "BM after feedback to cra is removed")+
  scale_color_manual(values = c("1" = "turquoise4", "0" = "firebrick2"))+
  scale_shape_manual(values=c(3, 16))+
  geom_text_repel(
    data = subset(Ecoli_Abla_df, Gene %in% labs_req),
    aes(label = Gene),
    vjust = -1, # adjust vertical position
    position = pos,
    color = "darkmagenta", # optional highlight color
    min.segment.length = 0,segment.color = NA
  )+stat_cor(method="spearman")


ggsave("Abla_cra.pdf")
ggsave("Abla_cra.jpeg")


####################################### 10 feedback metabolic genes

not_req_MGR_genes <- c("crp" , "fur" , "cra" , "lrp",  "pdhR", "argR" ,"purR" ,"cysB", "fhlA", "nagC")
aa <- match(not_req_MGR_genes,colnames(dist_matrix_1))

cc <- c()
for(i in 1:length(aa)){
  bb <- rownames(dist_matrix_1[dist_matrix_1[[aa[i]]] == 1,])
  cc <- c(cc,bb)
}



cc_ids <- match(cc,Ecoli_Abla_df$Gene)
cc_shp <- ifelse(1:nrow(Ecoli_Abla_df) %in% cc_ids, "1_hop", "Normal")

Ecoli_Abla_df$hop_1_10 <- cc_shp


labs_req <- Ecoli_Abla_df[Ecoli_Abla_df$CF_final_labels ==1 & Ecoli_Abla_df$Abla_10_BM_Final_labs ==0,1]  #"fdnI" "ddpF" "ddpD" "ddpC" "ddpB" "ddpA" "ddpX" "astC" "gabT" "cysA" "serC" "cysW" "cysP" "cysU" "argD"
highlight_ids <- match(labs_req, Ecoli_Abla_df$Gene)


intersect(labs_req,cc)  ##  "astC" "gabT" "cysA" "serC" "cysW" "cysP" "cysU" "argD"

pos <- position_jitter(width = 0.1, height = 0.1)

ggplot(Ecoli_Abla_df, aes(CF_S_pred, Abla_10_BM_Final))+geom_jitter(position = pos,aes(color = as.factor(True_labels), shape = as.factor(hop_1_10)), size = 2, alpha = 0.4)+geom_abline()+theme_bw()+
  labs(title = "Biomass of gene KOs before vs after feedback to top \n10 genes (highest 1-hop) are removed",x = "BM before feedback to top 10 genes (highest 1-hop) are removed", y = "BM after feedback to top 10 genes (highest 1-hop) are removed")+
  scale_color_manual(values = c("1" = "turquoise4", "0" = "firebrick2"))+
  scale_shape_manual(values=c(3, 16))+
  geom_text_repel(
    data = subset(Ecoli_Abla_df, Gene %in% labs_req),
    aes(label = Gene),
    vjust = -3, # adjust vertical position
    position = pos,
    color = "darkmagenta", # optional highlight color
    min.segment.length = 0,segment.color = NA
  )+stat_cor(method="spearman")


ggsave("Abla_10_metfeedbackgenes.pdf")
ggsave("Abla_10_metfeedbackgenes.jpeg")



####################################### 40 feedback metabolic genes

not_req_MGR_genes <- c("crp" , "fur" , "cra" , "lrp",  "pdhR", "argR" ,"purR" ,"cysB", "fhlA", "nagC",
                       "metJ" ,"paaX" ,"cytR", "argP", "cbl" , "glpR", "tyrR","malT", "gntR", "araC",
                       "trpR", "galR", "galS" ,"allR","exuR" ,"deoR", "mhpR" ,"lsrR" ,"hcaR", "idnR" ,
                       "uxuR" ,"rbsR", "xylR", "fucR","prpR" ,"caiF", "metR" ,"rhaS" ,"nanR" ,"rutR")

aa <- match(not_req_MGR_genes,colnames(dist_matrix_1))

cc <- c()
for(i in 1:length(aa)){
  bb <- rownames(dist_matrix_1[dist_matrix_1[[aa[i]]] == 1,])
  cc <- c(cc,bb)
}



cc_ids <- match(cc,Ecoli_Abla_df$Gene)
cc_shp <- ifelse(1:nrow(Ecoli_Abla_df) %in% cc_ids, "1_hop", "Normal")

Ecoli_Abla_df$hop_1_40 <- cc_shp


labs_req <- Ecoli_Abla_df[Ecoli_Abla_df$CF_final_labels ==1 & Ecoli_Abla_df$Abla_40_BM_Final_labs ==0,1]  #"astC" "gabT" "cysA" "serC" "cysW" "cysP" "cysU" "argD"
highlight_ids <- match(labs_req, Ecoli_Abla_df$Gene)


intersect(labs_req,cc)  ##  "astC" "gabT" "cysA" "serC" "cysW" "cysP" "cysU" "argD"

pos <- position_jitter(width = 0.1, height = 0.1)

ggplot(Ecoli_Abla_df, aes(CF_S_pred, Abla_MGR_40_BM_Final))+geom_jitter(position = pos,aes(color = as.factor(True_labels), shape = as.factor(hop_1_40)), size = 2, alpha = 0.4)+geom_abline()+theme_bw()+
  labs(title = "Biomass of gene KOs before vs after feedback to top \n40 genes (highest 1-hop) are removed",x = "BM before feedback to top 40 genes (highest 1-hop) are removed", y = "BM after feedback to top 40 genes (highest 1-hop) are removed")+
  scale_color_manual(values = c("1" = "turquoise4", "0" = "firebrick2"))+
  scale_shape_manual(values=c(3, 16))+
  geom_text_repel(
    data = subset(Ecoli_Abla_df, Gene %in% labs_req),
    aes(label = Gene),
    vjust = -1, # adjust vertical position
    position = pos,
    color = "darkmagenta", # optional highlight color
    min.segment.length = 0,segment.color = NA
  )+stat_cor(method="spearman")


ggsave("Abla_40_metfeedbackgenes.pdf")
ggsave("Abla_40_metfeedbackgenes.jpeg")

















####################

#install.packages("caret", type = "binary")
library(caret)


overall_precision_vec <- c()
overall_recall_vec <- c()
overall_F1_vec <- c()
overall_BA_vec <- c()




cx <- confusionMatrix(as.factor(Ecoli_Abla_df$CF_final_labels),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])



cx <- confusionMatrix(as.factor(Ecoli_Abla_df$Abla_crp_BM_Final_labs),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])


cx <- confusionMatrix(as.factor(Ecoli_Abla_df$Abla_fur_BM_Final_labs),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])


cx <- confusionMatrix(as.factor(Ecoli_Abla_df$Abla_cra_BM_Final_labs),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])


cx <- confusionMatrix(as.factor(Ecoli_Abla_df$Abla_10_BM_Final_labs),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])



cx <- confusionMatrix(as.factor(Ecoli_Abla_df$Abla_40_BM_Final_labs),as.factor(Ecoli_Abla_df$True_labels),positive = "0")

overall_precision_vec <- c(overall_precision_vec,cx$byClass["Precision"])
overall_recall_vec <- c(overall_recall_vec,cx$byClass["Recall"])
overall_F1_vec <- c(overall_F1_vec,cx$byClass["F1"])
overall_BA_vec <- c(overall_BA_vec,cx$byClass["Balanced Accuracy"])



abla_df <- data.frame(
  "Model_Type" = c("Model w/o any ablation","Model with crp ablation","Model with fur ablation","Model with cra ablation","Model with 10 metabolic feedback genes ablation","Model with 40 metabolic feedback genes ablation"),
  "Precision" = overall_precision_vec,
  "Recall" = overall_recall_vec,
  "F1_score" = overall_F1_vec,
  "Balanced Accuracy" = overall_BA_vec
)

write.csv(abla_df,"Metrics_for_Ablated_models.csv")


