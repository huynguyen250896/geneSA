# geneSA v0.1.1
#### I. Introduction
---
The package geneSA is built to serve as a support tool for the paper "*[Improving existing analysis pipeline to identify and analyze cancer driver genes using multi-omics data](https://www.nature.com/articles/s41598-020-77318-1)*". </br> A log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine an association between each gene and the survival rates of patients separately, and then adjust identified log-rank P-values following Benjamini-Hochberg FDR. Genes with adjusted log-rank P-values (also known as Q-values) <= 0.05 are preserved. </br> 

#### II. Understanding the tool
---
The following are parameters provided by geneSA: </br> 
data: frame or matrix. It represents its rows are genomic features and its columns are samples.
Note that samples in rows of `data` are included in your clinical data and in exactly the same order.
  
time: numeric or integer column vector. It is overall survival time of all samples extracted from 
your clinical data. Note that samples in rows of clinical data are included in `data` and in exactly the 
same order before extracting it.

status: binary column vector. It is overall survival status of all samples extracted from your clinical 
data (usually coded as 1 = death, 2 = alive). Note that samples in rows of clinical data are included in `data` and 
in exactly the same order before extracting it.

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
geneSA(data = exp1, time = clinical_exp$OS_MONTHS, status = clinical_exp$status)
```
#### V. Citation
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
Reference Type: Journal Article
Author: Nguyen, Quang-Huy
Le, Duc-Hau
Year: 2020
Title: Improving existing analysis pipeline to identify and analyze cancer driver genes using multi-omics data
Journal: Scientific Reports
Volume: 10
Issue: 1
Pages: 20521
Date: 2020/11/25
ISSN: 2045-2322
DOI: 10.1038/s41598-020-77318-1
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
