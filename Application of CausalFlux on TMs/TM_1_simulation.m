# ########################  Application for the same 
library(readxl)
library(matlabr)
library(writexl)

############################################  Functions for simulations

Simulate_TM_1 <- function(wd,e,c){
  
  setwd(wd)
  
  write.csv(e, file = "Exch_b.csv", row.names = FALSE, col.names = FALSE)
  write.csv(c, file = "Condition.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  matlabr::run_matlab_script("TM_1_simulation.m", display = TRUE, verbose = TRUE)
  
}

Simulate_TM_2 <- function(wd,e,c){
  
  setwd(wd)
  
  write.csv(e, file = "Exch_b.csv", row.names = FALSE, col.names = FALSE)
  write.csv(c, file = "Condition.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  matlabr::run_matlab_script("TM_2_simulation.m", display = TRUE, verbose = TRUE)
  
}

Simulate_TM_3 <- function(wd,e,c){
  
  setwd(wd)
  
  write.csv(e, file = "Exch_b.csv", row.names = FALSE, col.names = FALSE)
  write.csv(c, file = "Condition.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  matlabr::run_matlab_script("TM_3_simulation.m", display = TRUE, verbose = TRUE)
  
}


#############################################  calling the functions for simulations

wd <- "D:/work/Integrated_network_model/TMS_actual_simulations/"
Simulate_TM_2(wd,3.2,"WT")















