# Modeling and dissecting bidirectional feedback in gene-metabolite systems using the CausalFlux (CF) method

This is the official repository of the paper "Modeling and dissecting bidirectional feedback in gene-metabolite systems using the CausalFlux method" by Nilesh Subramanian, Pavan Kumar, Raghunathan Rengasamy, Nirav Bhatt, and Manikandan Narayanan.
 <!-- ![GH_1-1](https://github.com/BIRDSgroup/CausalFlux/blob/main/Application%20of%20CF%20on%20TMs/GH_1-1.png) -->
<img src="https://github.com/BIRDSgroup/CausalFlux/blob/main/Application%20of%20CF%20on%20TMs/GH_1-1.png" alt="CausalFlux Overview" style="width:35%; height:auto;">

## License Preamble
Copyright 2025 BIRDS Group, IIT Madras

CF framework is a free pipeline: you can redistribute it and modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

CF framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please take a look at the GNU Lesser General Public License for more details.

## This repository contains three sections:

1. CausalFlux (CF) framework
2. Application of CF framework on the testbed models (TMs)
3. Application of CF framework on the real-world data (E. coli)

### Section 1: Our CF framework
This section presents the code for the [CF framework](https://github.com/BIRDSgroup/CausalFlux/tree/main/Our%20CF%20framework), applicable to any generic case. It includes information on the datasets required and the function arguments. A detailed methodology is described in our paper.

### Section 2: Application of CF framework on the testbed models (TMs)
This section provides the following information:
1) Codes/data for simulating the steady-state fluxes and gene expression data in the three TMs for the WT/KO cases under all the exchange rates.
2) Codes/data for running the CF on the TMs [Code and Datasets](https://github.com/BIRDSgroup/CausalFlux/tree/main/Application%20of%20CF%20on%20TMs/Applying%20methods%20on%20TMs/CF/Code%20and%20Datasets).
3) Codes/data for generating the necessary figures in the manuscript and supplemenatry material[Code and Datasets](https://github.com/BIRDSgroup/CausalFlux/tree/main/Application%20of%20CF%20on%20TMs/Applying%20methods%20on%20TMs/CF/Code%20and%20Datasets).

### Section 3: Application of CF framework on the real-world data (E. coli)
This section provides the following information:
1) Codes/data for running the CF to perform single gene KOs in _E. coli_ [CausalFlux Runs](https://github.com/BIRDSgroup/CausalFlux/tree/main/Application%20of%20CF%20on%20E.%20coli/CausalFlux_Runs)
2) Codes/data to reconstruct/learn the parameters of the GRN (Gene Regularoty Network) for _E. coli_


## Glossary
- CF - Causal Flux (our methodology)
- TM - Toy Model
- GRN - Gene Regulatory Network
- KO - Knock out

## Session information 
### R 
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 26100)

R Packages required:
readxl, matlabr, bnlearn, dplyr, writexl, tidyverse , caret, ggplot2, ggpubr, ggExtra, gridExtra, ff
(ff_4.0.12, bit_4.0.5, gridExtra_2.3, ggExtra_0.10.1, ggpubr_0.6.0, caret_6.0-94, lattice_0.20-45, lubridate_1.9.2, forcats_1.0.0, stringr_1.5.0, purrr_1.0.1, readr_2.1.4,        
tidyr_1.3.0, tibble_3.2.1, ggplot2_3.4.4, tidyverse_2.0.0, writexl_1.4.2, dplyr_1.1.2, bnlearn_4.8.1, matlabr_1.5.2, readxl_1.4.2, graph_1.76.0, BiocGenerics_0.44.0)

### MATLAB
9.12.0.2039608 (R2022a) Update 5
COnstraint-Based Reconstruction and Analysis The COBRA Toolbox - 2026
Solver: Gurobi

