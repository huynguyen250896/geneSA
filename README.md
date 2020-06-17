# geneSA v0.1.0
#### I. Introduction
The package geneSA is built to serve as a support tool for the paper "An improved protocol for identification and analysis of driver genes using multi-omics data". A log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between the expression of each gene and the survival rates of patients, and then compute Q-value (Benjamini-Hochberg procedure) using the function computeQ for the package [computeQ](https://github.com/huynguyen250896/computeQ) based on the previously identified p-values. Genes with Q-value <= 0.05 are preserved. </br> 

#### II. Data Struture 
You must preprare the two kinds of the following data: *vt* and *df* (see the 'III.Implementation' section). </br> 
vt: a vector comprises genes of interest that you want to perform an survival analysis association with them individually. </br> 
df: a data frame comprises its rows are patients, its columns are genes of interest (the same as *vt*). NOTE that, the two last columns must essentially be survival times of patients and survival status of patients (binary variable). </br> 
Please download datasets [Dataset.zip](https://github.com/huynguyen250896/geneSA/blob/master/Dataset.zip) as examples to well grasp GeneSA's requirement on data structure. </br> 

#### III. Implementation
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/geneSA")
```
Call the library;
```sh
library(geneSA)
```
running example:
```sh
geneSA(genename = vt, event = df)
```
#### IV. Citation
Please kindly cite the following paper if you use the code, datasets or any results in this repo: </br>
```sh
...
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
