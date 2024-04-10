#' Calculation of differential expression for interacting genes of dysregulated edges between mutant and wild type samples 
#'
#' @param  exp.profile mRNA expression profile with row as genes, column as samples.
#' @param  Dysnet The dyregulated edges for each sample.
#' @param  Mut A gene mutation matrix with row as samples, column as genes. The values of "1" and "0"
#' representing mutant (case) and wild type (control) regarding the expression profile samples. 
#' @param  pvalueCutoff  The threshold of p value.   
#'
#' @return A list containing gene expression differences of dysregulated edges for each mutated gene

wilcoxRes <- function(exp.profile, Dysnet, Mut, pvalueCutoff = 0.05) {
  
  ###
  #' 0. Mutated sample labels ("1": case, "0": control)
  #'
  ################################
  sample.groups <- NULL
  Mutgenes <- NULL
  for(j in 1:dim(Mut)[2]){
    y1=which(Mut[,j]==1)
    if(length(y1)>0){
      sample.group <- rep("0",dim(Mut)[1])
      sample.group[y1] <- rep("1",length(y1))
      sample.groups <- cbind(sample.groups,sample.group)
      Mutgenes <- c(Mutgenes,colnames(Mut)[j])
    }
  }
  colnames(sample.groups) <-Mutgenes 
  rownames(sample.groups) <-rownames(Mut)
  
  
  ###
  #' 1. Gene expression profile of perturbed interactions
  #'
  ################################
  Dysgenes <- union(Dysnet$GeneA,Dysnet$GeneB)
  intergene <- intersect(Dysgenes,rownames(exp.profile))
  intersample <- intersect(rownames(sample.groups),colnames(exp.profile))
  
  DysExp <- exp.profile[intergene,intersample]
  DysExps <- 2^DysExp-1
  
  # For each mutated gene, calculate gene expression differences between MUT and WT samples
  wilcoxRes <- lapply(1:dim(sample.groups)[2],function(j) { 
    if(j%%10==0) {print(j)}
    sample.group <-sample.groups[,j] 
    
    ###
    #' 1. Calculate log2FC  of gene expression differences between MUT and WT samples
    #'
    ################################
    
    # Calculate the mean of sample expression
    group.mean <- tapply(1:ncol(DysExps), sample.group, function(sub.exp) {
      tmp.mean <- rowMeans(as.matrix(DysExps[, sub.exp]))
      return(tmp.mean)
    })
    # Avoid log2(0) or 0 in the denominator
    group.mean[["1"]][which(group.mean[["1"]] == 0)] <- 1e-05
    group.mean[["0"]][which(group.mean[["0"]] == 0)] <- 1e-05
    
    log2FC <- log2(group.mean[["1"]] / group.mean[["0"]])
    
    
    ###
    #' 2. Compare the expression of each dysgene between mutant and WT samples using Wilcoxon rank sum test
    #'
    ################################
    
    data.merge <- cbind.data.frame(t(DysExp), sample_group = sample.group)
    #replace gene name with specific symbol
    temp <- grep("-",colnames(data.merge))
    if(length(temp)>0){
      xx <- colnames(data.merge)[temp]
      colnames(data.merge)[temp] <- gsub("-","",xx)
      names(xx) <- gsub("-","",xx)
    }else{
      xx <- NULL
    }
    
    wilcox.res <- lapply(colnames(data.merge)[1:dim(DysExp)[1]], function(s) { 
      tmp.data <- data.merge[, c("sample_group", s)]
      tmp.data$sample_group <- factor(tmp.data$sample_group, levels = c("1", "0"))
      tmp.result <- coin::wilcox_test(as.formula(paste0(s, " ~ sample_group")), data = tmp.data)
      refine.res <- data.frame(signature = s, Zvalue = coin::statistic(tmp.result, type = "standardized"), p.value = coin::pvalue(tmp.result))
      return(refine.res)
    })
    wilcox.res <- do.call(rbind.data.frame, wilcox.res)
    colnames(wilcox.res)[2] <- "Zvalue"
    
    xxlocs <- match(names(xx),wilcox.res[,"signature"])
    wilcox.res[xxlocs,"signature"] <- xx
    
    sigwilcox <- wilcox.res[which(wilcox.res[,"p.value"] < pvalueCutoff),] 
    
    # return results
    res <- list(log2FC = log2FC, wilcox.res = wilcox.res, sigwilcox = sigwilcox,sample.group=sample.group)
    return(res)
  })
  
  names(wilcoxRes) <- colnames(sample.groups)
  
  # return results
  return(wilcoxRes)
  
}