setwd("D:/work/Integrated_network_model/Git_hub_codes/Ecoli/Application_other_media/")
aom <- read.csv("Ecoli_GR_pred_diff_media.csv", header = T, row.names = 1)


library(ggplot2)
ggplot(aom, aes(x = LB_CF_S, y = M9_CF_S)) +
  geom_jitter(
    color = "firebrick",     # Point color
    alpha = 0.4,             # Transparency
    size = 2,                # Point size
    width = 0.1, height = 0.1  # Jitter amount
  ) +
  geom_abline(linetype = "dashed", color = "red") +
  labs(
    title = "Predictions for LB media vs M9 media (jitter plot)",
    x = "Growth in LB media",
    y = "Growth in M9 media"
  ) +
  theme_bw()

ggsave("CF_lb_m9_preds.pdf")
ggsave("CF_lb_m9_preds.jpeg")



ggplot(aom, aes(x = LB_CF_S, y = TSB_CF_S)) +
  geom_jitter(
    color = "firebrick",     # Point color
    alpha = 0.4,             # Transparency
    size = 2,                # Point size
    width = 0.1, height = 0.1  # Jitter amount
  ) +
  geom_abline(linetype = "dashed", color = "red") +
  labs(
    title = "Predictions for LB media vs TSB media (jitter plot)",
    x = "Growth in LB media",
    y = "Growth in TSB media"
  ) +
  theme_bw()

ggsave("CF_lb_tsb_preds.pdf")
ggsave("CF_lb_tsb_preds.jpeg")



## zeros in LB

zero_LB <- aom[aom$LB_CF_S == 0,]  #63

zero_M9 <- aom[aom$M9_CF_S == 0,]  #133

zero_TSB <- aom[aom$TSB_CF_S == 0,]  #93


############################################################################# Bar plot using ggplot 

LB_gene <- zero_LB$Gene
M9_gene <- zero_M9$Gene
TSB_gene <- zero_TSB$Gene
common_gene <- intersect(intersect(LB_gene, M9_gene), TSB_gene)

# Create a data frame with the sizes
df <- data.frame(
  Set = c("LB media", "M9 media", "TSB media","Common genes"),
  Size = c(length(LB_gene), length(M9_gene), length(TSB_gene),length(common_gene))
)

# Bar plot using ggplot2
library(ggplot2)

ggplot(df, aes(x = Set, y = Size, fill = Set)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "CausalFlux predictions", y = "Number of predicted essential genes", x = " ")

ggsave("ess_gene_pred_CF.pdf")
ggsave("ess_gene_pred_CF.jpeg")



###########################################################################################
######## Checking for KOs in M9 and LB media

zeros_ind_CF_S_LB <- which(aom$LB_CF_S == 0) #63

zeros_ind_CF_S_M9 <- which(aom$M9_CF_S == 0) #133

zeros_ind_CF_S_TSB <- which(aom$TSB_CF_S == 0) #93

genes_CF_S_LB <- aom$Gene[zeros_ind_CF_S_LB]
genes_CF_S_M9 <- aom$Gene[zeros_ind_CF_S_M9]
genes_CF_S_TSB <- aom$Gene[zeros_ind_CF_S_TSB] 

inter_LB_M9_KO <- intersect(genes_CF_S_LB,genes_CF_S_M9) #63

inter_LB_TSB_KO <- intersect(genes_CF_S_LB,genes_CF_S_TSB) #61

inter_TSB_M9_KO <- intersect(genes_CF_S_TSB,genes_CF_S_M9) #99

differ_LB_M9_KO <- setdiff(genes_CF_S_M9,genes_CF_S_LB) #70

differ_LB_TSB_KO <- setdiff(genes_CF_S_TSB,genes_CF_S_LB) #32

differ_TSB_M9_KO <- setdiff(genes_CF_S_M9,genes_CF_S_TSB) #45





col1 <- c("LB essential","M9 essential", "TSB essential","Common","growth_in_LB_not_M9","growth_in_TSB_not_M9")
col2 <- c(df$Size[1],df$Size[2],df$Size[3],df$Size[4],length(differ_LB_M9_KO),length(differ_TSB_M9_KO))


col3 <- c(0,0,0,0,0,0)

df_suppl <- data.frame(col1,col2,col3)
colnames(df_suppl) <- c("Description","Number of genes","List of genes")

df_suppl$`List of genes`[[1]] <- paste(genes_CF_S_LB, collapse = ",")
df_suppl$`List of genes`[[2]] <- paste(genes_CF_S_M9, collapse = ",")
df_suppl$`List of genes`[[3]] <- paste(genes_CF_S_TSB, collapse = ",")
df_suppl$`List of genes`[[4]] <- paste(common_gene, collapse = ",")
df_suppl$`List of genes`[[5]] <- paste(differ_LB_M9_KO, collapse = ",")
df_suppl$`List of genes`[[6]] <- paste(differ_TSB_M9_KO, collapse = ",")

df_suppl <- as.data.frame(df_suppl)

write.csv(df_suppl, "ST3.csv")











