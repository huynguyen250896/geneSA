# geneSA v0.1.5
#### I. Introduction
---
The package geneSA is built to serve as a support tool for the paper "*[Improving existing analysis pipeline to identify and analyze cancer driver genes using multi-omics data](https://www.nature.com/articles/s41598-020-77318-1)*". </br> A log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine an association between each gene and the survival rates of patients separately, and then adjust identified log-rank P-values following Benjamini-Hochberg FDR. Genes with adjusted log-rank P-values (also known as Q-values) <= 0.05 are preserved. </br> 

#### II. Understanding the tool
---
The following are parameters provided by geneSA:

- data: data frame or matrix. It represents its rows are samples and its columns are genomic features .
Note that samples in rows of `data` are also included in your clinical data and in exactly the same order.
Users can feed any -omics data to this parameter; e.g., gene expression, copy number alteration,
methylation, or the like. NOTE that you must code values/observations in `data` as labels in advance. 
For example, expresison levels of genes usually are divided into highly expressed genes ("up") and lowly 
expressed genes ("down"), or into highly expressed genes ("up"), moderately expressed genes ("mid"), and 
lowly expressed genes ("down").
  
- time: numeric or integer column vector. It is overall survival time of all samples extracted from 
your clinical data. Note that samples in rows of clinical data are included in `data` and in exactly the 
same order before extracting it.

- status: binary column vector. It is overall survival status of all samples extracted from your clinical 
data (usually coded as 1 = death, 0 = alive). Note that samples in rows of clinical data are included in `data` and 
in exactly the same order before extracting it.

- Pcut: numeric. A user-defined P-value threshold to define statistical significance level. Default value is 
P-value <= 0.05

- Qcut: numeric. A user-defined Q-value threshold to define statistical significance level. Default value is 
Q-value <= 0.05

Please see & download data [data_n_code](https://github.com/huynguyen250896/geneSA/tree/master/data_n_code) as examples to well grasp the GeneSA's requirement
on data structure and its usage. </br> 

#### III. Pipeline
---
![Figure](https://imgur.com/hLlsaSl.png)
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
# exp is a matrix whose rows are samples and columns are genomic features
#>median is up-regulated genes and <median is down regulated genes
exp1 <- apply(exp,2, function(x) ifelse(x > median(x),"up","down")) %>% as.data.frame()

#Make sure samples that in rows of exp1 are also included in rows of clinical_exp and in exactly the same order
all(rownames(exp1) == rownames(clinical_exp))
#[1] FALSE
exp1 = exp1[rownames(clinical_exp),]

#RUN!!!
geneSA(data = exp1, time = clinical_exp$OS_MONTHS, status = clinical_exp$status, Pcut = 0.05, Qcut= 0.05)
```

#### V. What's new
---
- 2021-12-21: Now users have been able to define statistical significance level of their choice in relation to P-value and Q-value using `Pcut` and `Qcut`, respectively. 
- 2021-01-28: Now users have been able to divide the expression levels of genes over patients/samples into either two groups ("up" and "down") or three groups ("up", "mid", and "down")
- 2021-01-22: To be convenient more, I have changed the name of parameters in the geneSA. 

#### VI. Citation
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
