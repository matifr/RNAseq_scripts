##------------------------------------------
## Author: Matina Fragkogianni
## Date: 16-10-2019
##
## A function to make heatmaps
##
## @param matrix a matrix of gene expression 
## @param condition factorized vector 
## @param population factorized vector 
## @param cl_method clustering method to be used
## @ cl_distance_cols clustering metric for the columns
## @return heatmap 
##
##------------------------------------------

my_heatmap <- function(matrix, condition, population, cl_method = "complete", cl_distance_cols = "correlation", cl_cols = FALSE){
  
  require(pheatmap)
  df <- data.frame(Condition = condition, Population = population)
  rownames(df) <- colnames(matrix)
  names(df) <- c("Condition", "Population") 
  
  Condition        <- brewer.pal(3,"Dark2")
  names(Condition) <- levels(condition)
  
  Population        <- brewer.pal(4,"Set3")
  names(Population) <- levels(population)
  anno_colors <- list(Condition = Condition, Population = Population)
  
  pheatmap(matrix, annotation_col = df, annotation_colors = anno_colors,  scale = "row",
           show_rownames = TRUE, clustering_method = cl_method, clustering_distance_cols = cl_distance_cols, 
            show_colnames = TRUE, cluster_cols = cl_cols)
}


















