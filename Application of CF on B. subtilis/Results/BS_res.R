setwd("set it to directory where this csv file is stored")
tot_bs_df <- read.csv("Tot_BS_df.csv", header = T)


library(caret)

cxC <- caret::confusionMatrix(as.factor(tot_bs_df$CF_labs),as.factor(tot_bs_df$GT),positive = "0")

precision <- cxC$byClass["Precision"]
precision
recall <- cxC$byClass["Recall"]
recall
f1 <- 2 * (precision * recall) / (precision + recall)
f1
cxC
