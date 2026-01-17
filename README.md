# CausalFlux: an iterative causal surgery and linear optimization approach to model feedback in a gene-metabolite system

This is the official repository of the paper "CausalFlux: an iterative causal surgery and linear optimization approach to model feedback in a gene-metabolite system" by Nilesh Subramanian, Pavan Kumar, Raghunathan Rengasamy, Nirav Bhatt, and Manikandan Narayanan.
 <!-- ![GH_1-1](https://github.com/BIRDSgroup/CausalFlux/blob/main/Application%20of%20CF%20on%20TMs/GH_1-1.png) -->
<img src="https://github.com/BIRDSgroup/CausalFlux/blob/main/Application%20of%20CF%20on%20TMs/GH_1-1.png" alt="CausalFlux Overview" style="width:35%; height:auto;">

## License Preamble
Copyright 2025 BIRDS Group, IIT Madras

CF framework is a free pipeline: you can redistribute it and modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

CF framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please take a look at the GNU Lesser General Public License for more details.

## This repository contains three sections:

1. CF framework
2. Application of CF framework on the toy models (TMs)
3. Application of CF framework on the real-world data (E. coli)

### Section 1: Our CF framework
This section presents the code for the CF framework, applicable to any generic case. It includes information on the datasets required and the function arguments. A detailed methodology is described in our paper.

### Section 2: Application of CF framework on the toy models (TMs)
This section provides the following information:
1) Codes/data for simulating the steady-state fluxes and gene expression data in the three TMs for the WT/KO cases under all the exchange rates.
2) Codes/data for running the CF, TRIMER, and GIMME approach on the TMs.
3) Codes/data for generating the necessary figures in the manuscript and supplemenatry material.

### Section 3: Application of CF framework on the real-world data (E. coli)
This section provides the following information:
1) Codes/data for running the CF to perform single gene KOs in _E. coli_
2) Codes/data to reconstruct/learn the parameters of the GRN (Gene Regularoty Network) for _E. coli_
3) Codes/data for running the CF to perform double gene KOs in _E. coli_
4) Codes/data for generating the necessary table (Table 2 in the manuscript)

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

