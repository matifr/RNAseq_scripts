##------------------------------------------
## Author: Matina Fragkogianni
## Date: 16-10-2019
##
## A function that loads, selectes specific specific cell lines and filters genes based 
## a threshold or quantile gene expression
##
## @param toMatch a string to match and select appripriate cell line expression
## @param threshold a value for gene expression filtering
## @param prob a probabitly to filter based on 
## @return a dataframe with the filtered gene expression values
##
##------------------------------------------

getCCLE_data <- function(toMatch = "PANCREAS", threshold = 7, prob = NULL){  
  
  mat <- read.csv("~/Desktop/CCLE_RNAseq_genes_rpkm_20180929.csv", stringsAsFactors = FALSE)
  mat <- as_tibble(mat)

  if(!is.null(toMatch)){
    pMat <- mat %>% dplyr::select("Description", matches(toMatch))
  }else{
    pMat <- mat
  }

  Matmeans <- rowMeans(pMat[,2:ncol(pMat)])
  
  # Filtering for threshold > 7 (RPKM)
  cat(sprintf("Filtering for threshold >= %s", threshold))
  indexes = which(Matmeans >= threshold)
  Matmeans = data.frame("expression" = Matmeans[indexes])
  Matmeans$Description <- pMat$Description[indexes]
    
  if(!is.null(prob)){
    cat(sprintf("Filtering for probability >= %s\n", prob))
    q4 = quantile(Matmeans$expression, probs = prob)
    print(q4)
    Matmeans <- Matmeans[Matmeans$expression >= q4,]
  }else{
    print("I am returning the unfiltered matrix..")
  }
  
  return(Matmeans)
}



































