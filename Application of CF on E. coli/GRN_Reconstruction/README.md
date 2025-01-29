# GRN reconstruction

There are two folders: 

(1) **Code**: Contains all the neccessary code/data to reconstruct _E. coli_ GRN

(i) Load the "DAG_WL.csv" (that contains the TF-TG prior knowledge pairs); Load the "GE.rds" (contains the binarized gene expression data)

(ii) Modify the lines 1 and 2 in the "GRN_Reconstruction.R" script to load the "DAG_WL.csv" and "GE.rds" data correctly

(iii) Modify lines 7 and 12 in the "GRN_Reconstruction.R" script to save the output in the directory/folder of your choice

(2) **Results**: Contains the reconstructed GRN (bn_TR.RDS) and parameter estimation for the learnt GRN (bn_params_TR.RDS)
Load the "R" objects to view the results
