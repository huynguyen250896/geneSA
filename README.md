# geneSA v0.1.1
#### I. Introduction
---
The package geneSA is built to serve as a support tool for the paper "*Improving existing analysis pipeline to identify and analyze cancer driver genes using multi-omics data*". </br> A log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine an association between each gene and the survival rates of patients separately, and then adjust identified log-rank P-values following Benjamini-Hochberg FDR. Genes with adjusted log-rank P-values (also known as Q-values) <= 0.05 are preserved. </br> 

#### II. Data Struture 
---
You must preprare the two kinds of the following data: *vt* and *df* (see the 'III.Implementation' section). </br> 
vt: a vector comprises genes of interest that you want to perform a survival association analysis with them individually. </br> 
df: a data frame (e.g., gene expression data) comprises its rows are patients, its columns are genes of interest (the same as *vt*), its elements are gene expression levels (categorical variable: 1 = up-regulated, 0 = down-regulated, for example). NOTE that, the two last columns must essentially be (1) survival times of patients (continuous variables) and (2) survival status of patients (dichotomized variable: 1=death, 0 = alive). Alternatively, you can use DNA copy number alteration data or DNA methylation data or anything else (NOTE: categorical variable required) </br> 
Please download datasets [Dataset.zip](https://github.com/huynguyen250896/geneSA/blob/master/Dataset.zip) as examples to well grasp GeneSA's requirement on data structure. </br> 

#### III. Pipeline
---
![Figure](https://imgur.com/pvuJx9C.png)
**Figure:** Pipeline of the package geneSA.

#### IV. Implementation
---
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
#### V. Citation
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
...
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
