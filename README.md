# geneSA v0.1.0
#### I. Introduction
a log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between the expression of each gene and the survival rates of patients, and then compute Q-value (Benjamini-Hochberg procedure) using the function computeQ for the package [computeQ](https://github.com/huynguyen250896/computeQ) based on the previously identified p-values. Genes with Q-value <= 0.05 (Benjamini-Hocberg procedure) are preserved. </br> 


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
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3686388.svg)](https://doi.org/10.5281/zenodo.3686388)
```sh
@software{quang_huy_nguyen_2020_3686388,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/computeQ: v 0.1.0},
  month        = feb,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {0.1.0},
  doi          = {10.5281/zenodo.3686388},
  url          = {https://doi.org/10.5281/zenodo.3686388}
}
```
</br> And </br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3686391.svg)](https://doi.org/10.5281/zenodo.3829388)
```sh
@software{nguyen_quang_huy_2020_3829388,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/geneSA: GeneSA v0.1.0},
  month        = may,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.3829388},
  url          = {https://doi.org/10.5281/zenodo.3829388}
}
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
