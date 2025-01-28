# Intructions to run "CF_MTR_Ecoli.R" script

**Make sure all the folders and files in this folder are in the same working directory as that of the script**

## Modifications in "CF_MTR_Ecoli.R" script
1) Change the "curr_wd" (in line 1222) to the working directory where all the codes/files are from this folder. The results generated will be stored in this "curr_wd"
2) Change the directory in line 1225 to the directory where "**Ecoli_GT_keio_coolection.csv**" is present (this csv is currently in "/Application of CF on Ecoli/Results/Actual_GT").
3) Change the line in 977 to make sure the directory is pointed to the "GRN_RED_WO_BIGG" folder
4) Change the line in 984 to load the "GE.RData" in the "GRN_RED_WO_BIGG" folder
5) Change the line in 1015 to make sure the directory is pointed to the "MM_REQ" folder
   
## Modifications in the three ".m" scripts
1) In all three ".m" scripts, change line number 3 to the current working directory where all the codes/files from this folder are stored
2) Make sure the directory pointing to the "LB_media_constraints_iML1515.csv" (/MM/LB_media) is modified correctly in all the three ".m" scripts (line 106 in "CF_V2_Init_1_ID.m", "CF_V2_Iteration_2a_ID.m" and line 25 in " CF_V2_Iteration_2d_ID.m")


## For faster/parallel runs
1) Split the "**Ecoli_GT_keio_coolection.csv**" into multiple parts (say **n**)
2) Make **n** copies of this folder and change the necessary lines in the codes of each folder (see above) to reflect this change
3) Line 1225 from each of these **n** ".R" scripts should point to the **n** data splits respectively
4) Open **n** terminals and run the ".R" script from the corresponding folders
