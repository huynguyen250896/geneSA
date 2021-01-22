#' @title geneSA: Identification of genes significantly associated with patient outcome.
#'
#' @description a log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between each gene and patient outcome. Genes with Q-value <= 0.05 (Benjamini-Hocberg procedure) are preserved.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage computeC(data, time, status)
#'
#' @param data data frame or matrix. It represents its rows are genomic features and its columns are samples.
#' Note that samples in rows of \code{data} are included in your clinical data and in exactly the same order.
#'
#' @param time numeric or integer column vector. It is overall survival time of all samples extracted from 
#' your clinical data. Note that samples in rows of clinical data are included in \code{data} and in exactly the 
#' same order before extracting it.
#'
#' @param status binary column vector. It is overall survival status of all samples extracted from your clinical 
#' data (coded 1 = death, 2 = alive). Note that samples in rows of clinical data are included in \code{data} and 
#' in exactly the same order before extracting it.
#'
#' @return NULL
#'
#' @examples geneSA(data = exp1, time = clinical_exp$OS_MONTHS, status = clinical_exp$status)
#'
#' @export

geneSA = function(data = NULL, time = NULL, status = NULL){
  
  #Errors
  if(missing(data)){
    stop("Error: Input data is missing \n")
  }
  
  if(missing(time)){
    stop("Error: Overall survival time of all the patients is missing \n")
  }
  
  if(missing(status)){
    stop("Error: Overall survival status of all the patients is missing \n")
  }
  
  if(nrow(data) != length(time)){
    stop("Error: Please make sure samples in rows of data are included in your clinical data and in exactly the same order. \n")
  }
  
  if(nrow(data) != length(status)){
    stop("Error: Please make sure samples in rows of data are included in your clinical data and in exactly the same order. \n")
  }
  
  if(!(is.numeric(time) | is.integer(time))){
    stop("Error: Overall survival time must be numeric or integer. \n")
  }
  
  if(!(is.numeric(status) | is.integer(status))){
    stop("Error: Overall survival status must be numeric or integer. NOTE that status should be also as binary. \n")
  }
  
  #library
  library(survival)
  library(dplyr)
  
  #define the computeQ function, adjusting the log-rank P-values following Benjamini-Hochberg FDR
  computeQ <- function(x){
    (x$P.value*nrow(x))/(x$rank)
  }
  
  #run SA
  set.seed(25081996)
  dataset=cbind(data, time, status) %>% as.data.frame()
  df1=lapply(colnames(data),
             
             function(x) {
               
               formula <- as.formula(paste('Surv(time,status)~',as.factor(x)))
               coxFit <- survival::coxph(formula, data = dataset)
               summary(coxFit)
             })
  
  cc = data.frame(My_name_is = paste("Huy",1:length(df1)), HR=NA, confidence_intervals=NA, P.value=NA)
  
  for (i in c(1:length(df1))) {
    cc$HR[i] = round(df1[[i]][["coefficients"]][2],3) #hazard ratio
    cc$confidence_intervals[i] = paste(round(df1[[i]][["conf.int"]][[3]],3), "-", round(df1[[i]][["conf.int"]][[4]],3)) #95% CI
    cc$P.value[i] = df1[[i]][["logtest"]][3] #P-value
    rownames(cc)[i] =rownames(df1[[i]][["conf.int"]])
    order.pvalue = order(cc$P.value)
    cc = cc[order.pvalue,] #re-order rows following p-value
    cc$rank = c(1:length(df1)) #rank of P.value
    cc$Q.value = computeQ(cc) #compute Q-value
    rownames(cc) <- gsub("up","",rownames(cc)) #remove the word "up" in row names
  }
  cc=cc[,-1]
  cc = dplyr::select(cc, -rank) #remove the 'rank' column  
  cc = cc %>% subset(P.value <= 0.05) #only retain Genes with P <=0.05
  cc = cc %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
  write.table(cc,"gene_SA.txt",sep = "\t", quote = FALSE)
  
  #Messenge
  if(length(levels(as.factor(dataset[,1]))) <= 2){
    cat("\n","NOTE:" ,"\n","*gene_SA.txt placed in your current working directory.","\n","*Please check to identify which gene significantly associated with patient outcome.","\n","*In this case, the numerator is", levels(factor(dataset[,1]))[[2]], "and the denominator is", levels(factor(dataset[,1]))[[1]], ". In other words,", levels(factor(dataset[,1]))[[1]], "is considered as the reference group.")
  } else{
    for (j in 2:length(levels(as.factor(dataset[,1])))){
      cat("\n","NOTE:" ,"\n","*gene_SA.txt placed in your current working directory.","\n","*Please check to identify which gene significantly associated with patient outcome.","\n","*In this case, the numerator is", paste(levels(factor(dataset[,1]))[[i]], sep = ","), "and the denominator is", levels(factor(dataset[,1]))[[1]], ". In other words,", levels(factor(dataset[,1]))[[1]], "is considered as the reference group.")
    }}
}

