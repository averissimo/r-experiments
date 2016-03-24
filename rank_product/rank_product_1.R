### Rank product meta-analysis functions.
### To derive rank product p-values from a genes*samples data matrix 'dataMat', simply run rankProductWrapper(dataMat).


### Return a rank matrix in which [i, j] is the rank of object 'i' in condition 'j'.
### Assumes there are no NA's.
rankMatrix = function(mat)
{
  rankMat = mat
  for(i in 1:ncol(mat)){
    rankMat[sort(mat[, i], decreasing = T, index.return = T)$ix, i] = 1:nrow(mat)
  }
  return(rankMat)
}

### Sum of log-ranks.
### Takes as input a rank matrix.
logRankProduct = function(rankMat)
{
  return(apply(log(rankMat), 1, sum))
}

### Rank-product p-values computed via J. A. Koziol's approximation method.
### Takes as input an array of log-rank products (see logRankProduct function).
###
### References: J.A. Koziol. Comments on the Rank Product Method for Analyzing Replicated Experiments. FEBS Letters 584, 2010.
gammaPValues = function(lrp, nConds)
{
  nGenes = length(lrp)
  return(pgamma(-lrp + nConds * log(nGenes + 1), shape = nConds, lower.tail = F))
}

### Wrapper for rank product p-value computation. 
### Missing values are replaced by the minimum observed value.
### The input 'dataMat' is an expression matrix with dimensions #genes * #nSamples.
### Returns a list with the following fields:
### - rankMatrix: The rank matrix derived from dataMat. For every sample (column), each gene gets a rank from 1 to #genes.
### - logRankProduct: The log-rank product score derived for each gene.
### - pValues: The p-values resulting from the log-rank product scores.
rankProductWrapper = function(dataMat)
{
  dataMat[is.na(dataMat)] = min(dataMat)
  rankMat = rankMatrix(dataMat)
  lrp = logRankProduct(rankMat)
  pValues = gammaPValues(lrp, ncol(dataMat))
  return(list(rankMatrix = rankMat, logRankProduct = lrp, pValues = pValues))
}

######################################################################################

### Example with random matrix.
n = 5000 # Number of genes
m = 10 # Number of samples
numericMatrix = matrix(rnorm(n,m), n, m) # Create a random n  m matrix to simulate a differential expression data set.
output = rankProductWrapper(abs(numericMatrix)) # Compute the rank product p-values. NB: We sort genes according to their absolute differential expression.
outputRankMatrix = output$rankMatrix
outputLRP = output$logRankProduct
outputPValues = output$pValues