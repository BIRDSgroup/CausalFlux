library(caret)

Bin_Data_Ecoli <- data.frame(Act_vs_Pred_Ecoli$Gene,act_bin, pred_bin_cfs, pred_bin_cfmtr, pred_bin_tri)
colnames(Bin_Data_Ecoli) <- c("Genes","Binarized_Actual","Binarized_CFS_predictions","Binarized_CFMTR_predictions","Binarized_TRIMER_predictions")

write.csv(Bin_Data_Ecoli,"Binarized_act_pred_ecoli.csv", row.names = F)
