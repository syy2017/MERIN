#' Identify epigenetic regulators associated with immune gene sets
#'
#' @param exp.profile A numeric matrix containing the expression of mRNA with
#' rownames and colnames,  row as genes, column as samples.
#' @param Network The ligand-receptor interaction network. 
#' @param Mut A gene mutation matrix with row as samples, column as genes. The values of "1" and "0"
#' representing mutant (case) and wild type (control) regarding the expression profile samples.  
#' @param thr The significant level for detecting outliers 
#' @param pvalueCutoff  The threshold of p values to test differential expression. 
#' @param alpha The significance level used for identifying ER mutations related to dysregulated interactions.
#' @param n.sim Number of permutation times to calculate p values for ER mutations related to dysregulated interactions.
#'
#' @return A list including detailed correlation results for epigenetic regulator mutation-interaction pairs.

#' @examples
#' # test for "MERIPT" with example data
#' exp.profile <- load(file=paste0(workDir,"/data/PCG_tpm2.RData"))
#' Network <- load(file=paste0(workDir,"/Data/LRpairs.RData"))
#' Mut <- load(file=paste0(workDir,"/data/ER_driverMutMat.RData")) 
#' test_res <- MERIPT(exp.profile,Network,Mut)
#' test_res$MutPPI[1:3,] # showing the results
#' 
#' 


MERIPT <-  function(exp.profile, Network, Mut, thr=0.05, pvalueCutoff=0.05, alpha=0.05,n.sim=100){
  options (warn = -1)
  
  # Identify dysregulated ligand-receptor paris
  EdgeticDysRes <- EdgeticDys_CN(exp.profile,Network,thr)
  
  
  # Calculate expression difference
  Dysnet <- EdgeticDysRes[[1]] 
  wilcoxRess <- wilcoxRes(exp.profile, Dysnet, Mut, pvalueCutoff)
  
  
  
  #Identify epigenetic regulators associated with interaction dysregulation
  wilcoxRes <- wilcoxRess
  MutPurEdgeRes <- getEdgeticDrivers(wilcoxRes,Dysnet,Mut,alpha,n.sim)
  
  
  # return results
  return (MutPurEdgeRes)
}

