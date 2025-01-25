
library(Metrics)
library(caret)


Bin_Data_Ecoli <- data.frame(Act_vs_Pred_Ecoli$Gene,act_bin, pred_bin_cfs, pred_bin_cfmtr, pred_bin_tri)
colnames(Bin_Data_Ecoli) <- c("Genes","Binarized_Actual","Binarized_CFS_predictions","Binarized_CFMTR_predictions","Binarized_TRIMER_predictions")

write.csv(Bin_Data_Ecoli,"Binarized_act_pred_ecoli.csv", row.names = F)



x1 <- confusionMatrix(as.factor(Bin_Data_Ecoli$Binarized_CFS_predictions),as.factor(Bin_Data_Ecoli$Binarized_Actual))

balanced_accuracy_cfs <-  x1$byClass[11]

x2 <- confusionMatrix(as.factor(Bin_Data_Ecoli$Binarized_CFMTR_predictions),as.factor(Bin_Data_Ecoli$Binarized_Actual))

balanced_accuracy_cfmtr <-  x2$byClass[11]

x3 <- confusionMatrix(as.factor(Bin_Data_Ecoli$Binarized_TRIMER_predictions),as.factor(Bin_Data_Ecoli$Binarized_Actual))

balanced_accuracy_trimer <-  x3$byClass[11]

#### Recall

Recall_CF_s <- recall(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_CFS_predictions)) 
Recall_CF_mtr <- recall(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_CFMTR_predictions)) 
Recall_Trimer <- recall(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_TRIMER_predictions))

#### Precision

Precision_CF_s <- precision(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_CFS_predictions)) 
Precision_CF_mtr <- precision(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_CFMTR_predictions)) 
Precision_Trimer <- precision(as.factor(Bin_Data_Ecoli$Binarized_Actual),as.factor(Bin_Data_Ecoli$Binarized_TRIMER_predictions)) 

#### F1 score

f1_CF_s <- 2*(Recall_CF_s*Precision_CF_s)/(Recall_CF_s+Precision_CF_s)
f1_CF_mtr <- 2*(Recall_CF_mtr*Precision_CF_mtr)/(Recall_CF_mtr+Precision_CF_mtr)
f1_CF_Trimer <- 2*(Recall_Trimer*Precision_Trimer)/(Recall_Trimer+Precision_Trimer)

