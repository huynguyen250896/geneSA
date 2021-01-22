#library
if(!require(devtools)) install.packages("devtools")
devtools::install_github("huynguyen250896/geneSA")
library(geneSA)
library(tibble)
library(tidyverse)

#raw data
exp = read.table("data_mRNA_median_Zscores.txt", sep = '\t', check.names = F, header = T) 
clinical_exp = read.table("data_clinical_patient.txt", sep = '\t', check.names = F, header = T, fill = T)

#reduce dimension of exp to serve for application example
thirtygene=exp$Hugo_Symbol[19:30]

#only keep the 31 driver genes in exp and cna
exp=exp %>%
  dplyr::filter(.$Hugo_Symbol %in% thirtygene) %>%
  tibble::column_to_rownames('Hugo_Symbol') %>%
  dplyr::select(-Entrez_Gene_Id) %>% t()

#match patients that share between clinical data and EXP
clinical_exp = clinical_exp[clinical_exp$PATIENT_ID %in% rownames(exp),]

#create new column ‘status’ with survival event as binary
clinical_exp = clinical_exp %>%
  remove_rownames() %>%
  tibble::column_to_rownames('PATIENT_ID') %>%
  mutate(status = ifelse(clinical_exp$OS_STATUS == "LIVING",0,1)) #set event: death = 1, alive = 0

# create event vector for EXP
#>median is up-regulated genes and <median is down regulated genes
exp1 <- apply(exp,2, function(x) ifelse(x > median(x),"up","down")) %>% as.data.frame()
# check how many altered samples we have
table(as.matrix(exp1))
# down    up 
# 11425 11423 

#Make sure samples that in exp1 are also included in rows of clinical_exp and in exactly the same order
all(rownames(exp1) == rownames(clinical_exp))
#[1] FALSE
exp1 = exp1[rownames(clinical_exp),]

#What data look like
str(clinical_exp)
# 'data.frame':	1904 obs. of  21 variables:
#   $ LYMPH_NODES_EXAMINED_POSITIVE: int  10 0 1 3 8 0 1 1 1 0 ...
# $ NPI                          : num  6.04 4.02 4.03 4.05 6.08 ...
# $ CELLULARITY                  : chr  "" "High" "High" "Moderate" ...
# $ CHEMOTHERAPY                 : chr  "NO" "NO" "YES" "YES" ...
# $ COHORT                       : int  1 1 1 1 1 1 1 1 1 1 ...
# $ ER_IHC                       : chr  "Positve" "Positve" "Positve" "Positve" ...
# $ HER2_SNP6                    : chr  "NEUTRAL" "NEUTRAL" "NEUTRAL" "NEUTRAL" ...
# $ HORMONE_THERAPY              : chr  "YES" "YES" "YES" "YES" ...
# $ INFERRED_MENOPAUSAL_STATE    : chr  "Post" "Pre" "Pre" "Pre" ...
# $ INTCLUST                     : chr  "4ER+" "4ER+" "3" "9" ...
# $ AGE_AT_DIAGNOSIS             : num  75.7 43.2 48.9 47.7 77 ...
# $ OS_MONTHS                    : num  140.5 84.6 163.7 164.9 41.4 ...
# $ OS_STATUS                    : chr  "LIVING" "LIVING" "DECEASED" "LIVING" ...
# $ CLAUDIN_SUBTYPE              : chr  "claudin-low" "LumA" "LumB" "LumB" ...
# $ THREEGENE                    : chr  "ER-/HER2-" "ER+/HER2- High Prolif" "" "" ...
# $ VITAL_STATUS                 : chr  "Living" "Living" "Died of Disease" "Living" ...
# $ LATERALITY                   : chr  "Right" "Right" "Right" "Right" ...
# $ RADIO_THERAPY                : chr  "YES" "YES" "NO" "YES" ...
# $ HISTOLOGICAL_SUBTYPE         : chr  "Ductal/NST" "Ductal/NST" "Ductal/NST" "Mixed" ...
# $ BREAST_SURGERY               : chr  "MASTECTOMY" "BREAST CONSERVING" "MASTECTOMY" "MASTECTOMY" ...
# $ status                       : num  0 0 1 0 1 1 0 1 1 1 ...

dim(exp1)
#[1] 1904   12

#RUN!!!
geneSA(data = exp1, time = clinical_exp$OS_MONTHS, status = clinical_exp$status)


