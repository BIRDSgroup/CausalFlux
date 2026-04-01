setwd("D:/work/Integrated_network_model/Git_hub_codes/Ecoli/CausalFlux_TRIMER/")
oep <- read.csv("Overall_Ecoli_preds.csv", header = T, row.names = 1)


otp <- c()
ofp <- c()
ofn <- c()
otn <- c()
op <- c()
or <- c()
of <- c()
oba <- c()


library(caret)
cx <- confusionMatrix(as.factor(oep$CF_0th_labels),as.factor(oep$True_labels),positive = "0")


otp <- c(otp,cx$table[1,1])
ofp <- c(ofp,cx$table[1,2])
ofn <- c(ofn,cx$table[2,1])
otn <- c(otn,cx$table[2,2])
op <- c(op,cx$byClass["Precision"])
or <- c(or, cx$byClass["Recall"])
of <- c(of,  cx$byClass["F1"])
oba <- c(oba,  cx$byClass["Balanced Accuracy"])


cx <- confusionMatrix(as.factor(oep$CF_1st_labels),as.factor(oep$True_labels),positive = "0")

otp <- c(otp,cx$table[1,1])
ofp <- c(ofp,cx$table[1,2])
ofn <- c(ofn,cx$table[2,1])
otn <- c(otn,cx$table[2,2])
op <- c(op,cx$byClass["Precision"])
or <- c(or, cx$byClass["Recall"])
of <- c(of,  cx$byClass["F1"])
oba <- c(oba,  cx$byClass["Balanced Accuracy"])



cx <- confusionMatrix(as.factor(oep$CF_final_labels),as.factor(oep$True_labels),positive = "0")

otp <- c(otp,cx$table[1,1])
ofp <- c(ofp,cx$table[1,2])
ofn <- c(ofn,cx$table[2,1])
otn <- c(otn,cx$table[2,2])
op <- c(op,cx$byClass["Precision"])
or <- c(or, cx$byClass["Recall"])
of <- c(of,  cx$byClass["F1"])
oba <- c(oba,  cx$byClass["Balanced Accuracy"])


cx <- confusionMatrix(as.factor(oep$TR_labels),as.factor(oep$True_labels),positive = "0")

otp <- c(otp,cx$table[1,1])
ofp <- c(ofp,cx$table[1,2])
ofn <- c(ofn,cx$table[2,1])
otn <- c(otn,cx$table[2,2])
op <- c(op,cx$byClass["Precision"])
or <- c(or, cx$byClass["Recall"])
of <- c(of,  cx$byClass["F1"])
oba <- c(oba,  cx$byClass["Balanced Accuracy"])


somedf7 <- data.frame(
  "Methods" = c("CausalFlux (0th iteration)","CausalFlux (1st iteration)","CausalFlux (final iteration)","TRIMER"),
  "TP" = otp,
  "FP" = ofp,
  "FN" = ofn,
  "TN" = otn,
  "Precision" = op,
  "Recall" = or,
  "F1" = of,
  "Balanced_Accuracy" = oba
)

write.csv(somedf7,"Metrics_CausalFlux.csv")

