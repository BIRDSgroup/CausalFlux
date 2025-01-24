
setwd("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/MM_REQ/")
# Our IAF1260 data
met_gene_reg_data <- readRDS("Metabolite_Gene_Regulation.RDS")

Bigg_ids_gene_names <- readRDS("D:/work/Integrated_network_model/Ecoli_intg_ntwk/metabolic_aspect/Auto_RUN/MM_REQ/Bigg_ids_gene_names.RDS")
xids <- match(met_gene_reg_data$`TF-gene`, Bigg_ids_gene_names$Gene_names)


x_bigg <- Bigg_ids_gene_names$Bigg_ids[xids]

mgr <- met_gene_reg_data

mgr$`TF-gene` <- x_bigg

saveRDS(mgr, "MGR_bigg_ids.RDS")



























