#' @title Tool for Identification of genes significantly associated with prognostic value (survival rates of patients)
#'
#' @description a log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between the expression of each gene and the survival rates of patients. Genes with Q-value <= 0.05 (Benjamini-Hocberg procedure) are preserved.
#'
#' @param genename,event
#'
#' @return NULL
#'
#' @examples geneSA(cna = vt, exp1 = df)
#'
#' @export

geneSA = function(genename=NULL, event=NULL){

  #missing inputs
  if(missing(genename)){
    stop("Error: vector of genes considered is missing \n")
  }

  if(missing(event)){
    stop("Error: data frame of expression levels of genes with corresponding events is missing \n")
  }

  #library
  devtools::install_github("huynguyen250896/computeQ")
  if(!require(survival)) install.packages('survival')
  if(!require(rlist)) install.packages("https://cran.r-project.org/src/contrib/Archive/rlist/rlist_0.4.tar.gz", repos = NULL)
  library(computeQ)
  library('survival')
  library(rlist)
  library(dplyr)
  library(tidyverse)
  library(tidyr)

  #run SA
  set.seed(420)
  df1=lapply(genename,

             function(x) {

               formula <- as.formula(paste('Surv(time,event)~',as.factor(x)))
               coxFit <- coxph(formula, data = event)
               summary(coxFit)
             })


  df2 =  list.filter(df1, logtest[["pvalue"]] < 0.05) #preserve genes with P<=0.05

  cc = data.frame(No. = paste("Gene ",1:length(df2)), HR=NA, confidence_intervals=NA, P.value=NA)

  for (i in c(1:length(df2))) {
    cc$HR[i] = round(df2[[i]][["coefficients"]][2],3) #hazard ratio
    cc$confidence_intervals[i] = paste(round(df2[[i]][["conf.int"]][[3]],3), "-", round(df2[[i]][["conf.int"]][[4]],3)) #95% CI
    cc$P.value[i] = df2[[i]][["logtest"]][3] #P-value
    rownames(cc)[i] =rownames(df2[[i]][["conf.int"]])
    order.pvalue = order(cc$P.value)
    cc = cc[order.pvalue,] #re-order rows following p-value
    cc$rank = c(1:length(df2)) #rank of P.value
    cc$Q.value = computeQ(cc) #compute Q-value
    rownames(cc) <- gsub("up","",rownames(cc)) #remove the word "up" in row names
    write.table(cc,"gene_SA.txt",sep = "\t", quote = FALSE)
  }
  cc = read.table("gene_SA.txt", sep="\t", check.names = FALSE, row.names = 1, header = TRUE)
  cc = cc[,-c(1,5)] #remove the two unnescessary columns: No. and rank
  cc = cc %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
  View(cc)
  write.table(cc,"gene_SA.txt",sep = "\t", quote = FALSE)
  writeLines("NOTE:\n*gene_SA.txt placed in your current working directory\n*Please check to identify which genes significantly associated with prognostic value (survival rates of patients)\n*In any case, the numerator is up-expression level and the denominator is down-expression level. In other words, the down-expression level of genes considered is the reference.")
}
