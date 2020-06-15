# geneSA v0.1.0
#### I. Introduction
The package geneSA is built to serve as a support tool for the paper "An improved protocol for identification and analysis of driver genes using multi-omics data". A log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between the expression of each gene and the survival rates of patients, and then compute Q-value (Benjamini-Hochberg procedure) using the function computeQ for the package [computeQ](https://github.com/huynguyen250896/computeQ) based on the previously identified p-values. Genes with Q-value <= 0.05 are preserved. </br> 


#### IV. Implementation
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
Please kindly cite the two repositories if you use the code, datasets or any results in this repo: </br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3872669.svg)](https://doi.org/10.5281/zenodo.3872669)
```sh
@software{nguyen_quang_huy_2020_3872669,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/GeneCluster: GeneCluster v0.1.0},
  month        = jun,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.3872669},
  url          = {https://doi.org/10.5281/zenodo.3872669}
}
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
