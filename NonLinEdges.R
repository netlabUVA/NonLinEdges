###This function is the basic framework of the method presented in "Estimating nonlinear relations 
###among random variables: A network analytic approach." It differs from what is presented in that
###instead of residualizing on only the significantly linear edges, we residualize on regressions of
###all variables on all other variables. dat is the user's data in wide format, nperm is the number
###of permutations for the permutation test. The article uses 1000. The function returns the p-values
###of the partial distance correlations

NonLinEdges <- function(dat, nperm){
  ###Set up residuals dataframe
  residuals <- data.frame(matrix(ncol = ncol(dat), nrow = 0))
  colnames(resid) <- colnames(dat)
  ###Get residuals
  for(i in 1:ncol(residuals)){
    outcol <- colnames(residuals)[i]
    residuals[,i] <- resid(lm(outcol ~ .))
  }
  ###Get pdcor matrix
  k <- 1:ncol(residuals)
  p_mat <- matrix(0, nrow = ncol(residuals), ncol = ncol(residuals))
  return_cor <- data.frame(NoCN(residuals))
  for(i in 1:ncol(residuals)){
    for(j in 1:ncol(residuals)){
      k <- -c(i,j)
      test <- pdcor.test(residuals[,i], residuals[,j], residuals[,k], nperm)
      p_mat[i,j] <- test$p.value
    }
  }
  colnames(p_mat) <- colnames(residuals)
  rownames(p_mat) <- colnames(residuals)
  
  return(p_mat)

}