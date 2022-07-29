#' @title geneSA: Identification of genes significantly associated with patient outcome.
#'
#' @description a log-rank test in univariate Cox regression analysis with a proportional hazards model is performed to examine the association between each gene and patient outcome. Genes with Q-value <= 0.05 (Benjamini-Hocberg procedure) are preserved.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage geneSA(data, time, status, Pcut, Qcut)
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
#' @param Pcut numeric. A user-defined P-value threshold to define significance. default value is P-value <= 0.05.
#'
#' @param Qcut numeric. A user-defined Q-value threshold to define significance. default value is Q-value <= 0.05.
#'
#' @param univariate boolean. Whether geneSA runs an univariate or a multivariate survival analysis. Default value is univariate = T.
#'
#' @return NULL
#'
#' @examples geneSA(data = exp1, time = clinical_exp$OS_MONTHS, status = clinical_exp$status, univariate = T)
#'
#' @export

geneSA = function(data = NULL, time = NULL, status = NULL, Pcut = 0.05, Qcut = 0.05, univariate = T){
  
  #Error messages
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
    stop("Error: Overall survival status must be numeric or integer. \n")
  }
  
  #library
  library(survival)
  library(dplyr)
  
  #define the computeQ function, adjusting the log-rank P-values following Benjamini-Hochberg FDR
  computeQ <- function(x){
    (x$P.value*nrow(x))/(x$rank)
  }
  
  #Main function
  if(univariate == T | univariate == TRUE){ # univariate survival analysis
    if(length(table(as.matrix(data))) == 2){
      #run SA
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
        rownames(cc) <- gsub(levels(as.factor(dataset[,1]))[[2]],"",rownames(cc)) #remove the word "up" in row names
      }
      cc=cc[,-1]
      cc = dplyr::select(cc, -rank) #remove the 'rank' column  
      cc = cc %>% subset(P.value <= Pcut) #only retain Genes with P <=0.05
      cc = cc %>% subset(Q.value <= Qcut) #only retain Genes with Q <=0.05
      
      #Print
      write.table(cc,"gene_SA_univariate_2.txt",sep = "\t", quote = FALSE)
      
      #Messenger
      cat("\n","NOTE:" ,"\n","*gene_SA_univariate_2.txt stored in your current working directory.","\n","*Please check to identify which gene is significantly associated with patient outcome.","\n","*In this case, the numerator is", levels(factor(dataset[,1]))[[2]], "expression", "and the denominator is", levels(factor(dataset[,1]))[[1]], "expression.", " In other words,", levels(factor(dataset[,1]))[[1]], "expression", "is considered as the reference group.")
    } else{
      
      #run SA
      dataset=cbind(data, time, status) %>% as.data.frame()
      
      #extract value of k at which k^th column that has 3 levels assigned to value, otherwise assigned to value1
      value = data.frame(My_name_is = paste("Huy",1:ncol(data)), value = NA)
      value1 = data.frame(My_name_is = paste("Huy",1:ncol(data)), value1 = NA)
      for (k in 1:ncol(data)){
        if(length(table(data[,k])) == 3){
          value$value[[k]] <- k
        }else{
          value1$value1[[k]] <- k
        }
      }; value = na.omit(value); value1 = na.omit(value1)
      
      #relevel
      for (m in value[,2]){
        if(all(levels(factor(dataset[, value1[1,2]]))[[2]] != levels(factor(dataset[,m]))[[2]])){
          factor = c(levels(factor(dataset[,m]))[[1]], levels(factor(dataset[,m]))[[3]], levels(factor(dataset[,m]))[[2]])
          dataset[,m] <- factor(dataset[,m],
                                levels = factor)
        }
      }
      
      df1=lapply(colnames(data),
                 
                 function(x) {
                   
                   formula <- as.formula(paste('Surv(time,status)~',as.factor(x)))
                   coxFit <- survival::coxph(formula, data = dataset)
                   summary(coxFit)
                 })
      
      cc = data.frame(My_name_is = paste("Huy",1:length(df1)), HR1=NA, confidence_intervals1=NA, HR2=NA, confidence_intervals2=NA, P.value=NA)
      
      for (i in c(1:length(df1))) {
        if(length(levels(as.factor(dataset[,i]))) == 2){
          cc$HR1[i] = round(df1[[i]][["coefficients"]][2],3) #hazard ratio
          cc$confidence_intervals1[i] = paste(round(df1[[i]][["conf.int"]][[3]],3), "-", round(df1[[i]][["conf.int"]][[4]],3)) #95% CI
          cc$P.value[i] = df1[[i]][["logtest"]][3] #P-value
          rownames(cc)[i] =rownames(df1[[i]][["conf.int"]])
          order.pvalue = order(cc$P.value)
          cc = cc[order.pvalue,] #re-order rows following p-value
          cc$rank = c(1:length(df1)) #rank of P.value
          cc$Q.value = computeQ(cc) #compute Q-value
          rownames(cc) <- gsub(levels(as.factor(dataset[,i]))[[2]],"",rownames(cc)) #remove the word "up" in row names
          
        } else{
          cc$HR1[i] = round(df1[[i]][["coefficients"]][1,2],3) #hazard ratio 1
          cc$HR2[i] = round(df1[[i]][["coefficients"]][2,2],3) #hazard ratio 2
          cc$confidence_intervals1[i] = paste(round(df1[[i]][["conf.int"]][[1,3]],3), "-", round(df1[[i]][["conf.int"]][[1,4]],3)) #95% CI 
          cc$confidence_intervals2[i] = paste(round(df1[[i]][["conf.int"]][[2,3]],3), "-", round(df1[[i]][["conf.int"]][[2,4]],3)) #95% CI 
          cc$P.value[i] = df1[[i]][["logtest"]][3] #P-value
          rownames(cc)[i] =rownames(df1[[i]][["conf.int"]])[[2]]
          order.pvalue = order(cc$P.value)
          cc = cc[order.pvalue,] #re-order rows following p-value
          cc$rank = c(1:length(df1)) #rank of P.value
          cc$Q.value = computeQ(cc) #compute Q-value
          rownames(cc) <- gsub(levels(as.factor(dataset[,i]))[[2]],"",rownames(cc)) #remove the word "up" in row names
        }
      }
      
      cc=cc[,-1]
      cc = dplyr::select(cc, -rank) #remove the 'rank' column  
      cc = cc %>% subset(P.value <= Pcut) #only retain Genes with P <=0.05
      cc = cc %>% subset(Q.value <= Qcut) #only retain Genes with Q <=0.05
      
      #rename columns
      for (h in 1:ncol(data)){
        if(length(table(data[,h])) == 3){
          colnames(cc)[1] = paste0("HR_",levels(as.factor(dataset[,h]))[[2]])
          colnames(cc)[2] = paste0("confidence_intervals_",levels(as.factor(dataset[,h]))[[2]])
          colnames(cc)[3] = paste0("HR_",levels(as.factor(dataset[,h]))[[3]])
          colnames(cc)[4] = paste0("confidence_intervals_",levels(as.factor(dataset[,h]))[[3]])
        }
      }
      
      #Print 
      write.table(cc,"gene_SA_univariate_3.txt",sep = "\t", quote = FALSE)
      
      #Messenger
      cat("\\n","NOTE:" ,"\\n","*gene_SA_univariate_3.txt placed in your current working directory.","\\n","*Please check to identify which gene is significantly associated with patient outcome.","\\n","*In this case, the numerator is", paste(levels(factor(dataset[,value[1,2]]))[[2]], "and", levels(factor(dataset[,value[1,2]]))[[3]]), ", whereas the denominator is", levels(factor(dataset[,value1[1,2]]))[[1]],". In other words,", levels(factor(dataset[,value1[1,2]]))[[1]], "is considered as the reference group.")
    }
  } else{ # multivariate survival analysis
    #run SA
    dataset=cbind(data, time, status) %>% as.data.frame()
    
    
    coxFit <- survival::coxph(Surv(time,status)~ . , data = dataset)
    df1 = summary(coxFit)
    
    cc = data.frame(My_name_is = paste("Huy",1:nrow(df1[["coefficients"]])), HR=NA, confidence_intervals=NA, P.value=NA)
    
    for (i in c(1:nrow(df1[["coefficients"]]))) {
      cc$HR[i] = round(df1[["coefficients"]][i,2],3)
      cc$confidence_intervals[i] = paste(round(df1[["conf.int"]][i,3],3), "-", round(df1[["conf.int"]][i,3],3)) #95% CI
      cc$P.value[i] = df1[["coefficients"]][i,5] #P-value
      rownames(cc)[i] = rownames(df1[["coefficients"]])[i]
      order.pvalue = order(cc$P.value)
      cc = cc[order.pvalue,] #re-order rows following p-value
      cc$rank = c(1:nrow(df1[["coefficients"]])) #rank of P.value
      cc$Q.value = computeQ(cc) #compute Q-value
      rownames(cc) <- gsub(levels(as.factor(dataset[,1]))[[2]],"",rownames(cc)) #remove the word "up" in row names
    }
    
    cc=cc[,-1]
    cc = dplyr::select(cc, -rank) #remove the 'rank' column  
    cc = cc %>% subset(P.value <= Pcut) #only retain Genes with P <=0.05
    cc = cc %>% subset(Q.value <= Qcut) #only retain Genes with Q <=0.05
    
    #Print 
    write.table(cc,"gene_SA_multivariate.txt",sep = "\t", quote = FALSE)
    
    #Messenger
    cat("\n","NOTE:" ,"\n","*gene_SA_multivariate.txt stored in your current working directory.","\n","*Please check to identify which gene is significantly associated with patient outcome.","\n","*In this case, the numerator is", levels(factor(dataset[,1]))[[2]], "expression", "and the denominator is", levels(factor(dataset[,1]))[[1]], "expression.", " In other words,", levels(factor(dataset[,1]))[[1]], "expression", "is considered as the reference group.")
  }
}

