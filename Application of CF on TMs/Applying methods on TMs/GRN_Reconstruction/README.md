# Instructions for this section

This section provides instructions to run "GRN_TMS.R" script to reconstruct and learn the parameters of the GRN for the toy models.

"str_param()" function will perform GRN reconstruction and parameter learning.

## Instructions for "GRN_TMS.R" script

### TM1
1) Line 19 provides the whitelist edges for TM1 GRN (prior knowledge)

2) Lines from 21 to 29 should be run based on the genen expression data for the corresponding exchange rate.

   a) If the exchange rate is 3.2: Lines 22 and 23 must be excecuted. Remaining lines can be commented out. Change line 22 to working directory where the gene expression data corresponding to 3.2 exchange rate is stored.

   b) If the exchange rate is 320: Lines 25 and 26 must be excecuted. Remaining lines can be commented out. Change line 25 to working directory where the gene expression data corresponding to 320 exchange rate is stored.

   c) If the exchange rate is 3200: Lines 28 and 29 must be excecuted. Remaining lines can be commented out. Change line 28 to working directory where the gene expression data corresponding to 3200 exchange rate is stored.

3) Change line 34 to working directory to store the GRN as RDS 

### TM2
1) Line 41 provides the whitelist edges for TM2 GRN (prior knowledge)

2) Lines from 43 to 51 should be run based on the genen expression data for the corresponding exchange rate.

   a) If the exchange rate is 3.2: Lines 44 and 45 must be excecuted. Remaining lines can be commented out. Change line 44 to working directory where the gene expression data corresponding to 3.2 exchange rate is stored.

   b) If the exchange rate is 320: Lines 47 and 48 must be excecuted. Remaining lines can be commented out. Change line 47 to working directory where the gene expression data corresponding to 320 exchange rate is stored.

   c) If the exchange rate is 3200: Lines 50 and 51 must be excecuted. Remaining lines can be commented out. Change line 50 to working directory where the gene expression data corresponding to 3200 exchange rate is stored.

3) Change line 56 to working directory to store the GRN as RDS 

### TM3
1) Line 63 provides the whitelist edges for TM3 GRN (prior knowledge)

2) Lines from 65 to 73 should be run based on the genen expression data for the corresponding exchange rate.

   a) If the exchange rate is 3.2: Lines 66 and 67 must be excecuted. Remaining lines can be commented out. Change line 66 to working directory where the gene expression data corresponding to 3.2 exchange rate is stored.

   b) If the exchange rate is 320: Lines 69 and 70 must be excecuted. Remaining lines can be commented out. Change line 69 to working directory where the gene expression data corresponding to 320 exchange rate is stored.

   c) If the exchange rate is 3200: Lines 72 and 73 must be excecuted. Remaining lines can be commented out. Change line 72 to working directory where the gene expression data corresponding to 3200 exchange rate is stored.

3) Change line 77 to working directory to store the GRN as RDS


## Results from GRN reconstruction

The results (GRN reconstruction and parameters leanrt) for the toy models are in the TM1, TM2, and TM3 folders of "Applying methods on TMs/CF/CF_S/" in case of CF-S and "Applying methods on TMs/CF/CF_MTR/" in case of CF-MTR.




