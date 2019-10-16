calculate_correlation <- function(replicate1, replicate2){
  #caclulates correlation and prints it on the screen 
  cor = cor.test(as.matrix(replicate1), as.matrix(replicate2), method = "pearson")
  # plots the replicates on a scaterplot
  plot(as.matrix(replicate1), as.matrix(replicate2), xlab = "replicate1", ylab = "replicate2")
  return(cor)
}