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
2. Application of CF framework on the toy models
3. Application of CF framework on the real-world data

### Section 1: Our CF framework
This section presents the code for the CF-S and CF-MTR frameworks, applicable to any generic case. It includes information on the datasets required and the function arguments. A detailed methodology is described in our paper.

### Section 2: Application of CF framework on the toy models
This section provides the following information:
1) Codes/data for simulating the steady-state fluxes and gene expression data in the three TMs for the WT/KO cases under all the exchange rates.
2) Codes/data for running the CF-S/MTR, TRIMER, and GIMME approach on the TMs.
3) Codes/data for generating the necessary figures in the manuscript and supplemenatry material.

### Section 3: Application of CF framework on the real-world data
This section provides the following information:
1) Codes/data for running the CF-S/MTR tp perform single gene KOs in _E. coli_
2) Codes/data to reconstruct/learn the parameters of the GRN for _E. coli_
3) Codes/data for running the CF-S/MTR tp perform double gene KOs in _E. coli_
4) Codes/data for generating the necessary table (Table 2 in the manuscript)




