#' Detection of epigenetic regulator mutations involved in the dysregulated network edges
#
#'@param wilcoxRes Gene expression differences of dysregulated edges for each mutated gene.
#'@param Dysnet Dysregulated network for samples.
#'@param Mut A gene mutation matrix with row as samples, column as genes. The values of "1" and "0"
#' representing mutant (case) and wild type (control) regarding the expression profile samples.
#'@param alpha The significance level.
#'@param n.sim The simulation times to calculate p values.
#'
#'@return A list containing the mutation perturbed PPIs in cancer samples (MutPPI)


getEdgeticDrivers <- function(wilcoxRes,Dysnet,Mut,alpha=0.05,n.sim){
  
  #########################Mutation, PPI, #sample-1, #sample-2, #overlap, p
  
  N11=length(unique(Dysnet$Sam.ID))
  N22 <- 0
  for (i in 1:dim(Mut)[1]){
    xx <- Mut[i,]
    if(sum(xx)>0){
      N22 <- N22+1
    }
  }
  
  FF.nets <- lapply(1:length(wilcoxRes),function(j){
    if(j%%20==0) {print(j)}
    wilcoxRess <- wilcoxRes[[j]]
    sigwilcox <- wilcoxRess$sigwilcox
    sample.group <- wilcoxRess$sample.group
    mutgene <- names(wilcoxRes)[j] 
    
    #Mutated samples of the jth ER gene
    y1=which(sample.group=="1")
    if(length(y1)< 1){
      #
      return(NULL)
    }
    Mut.sam=unique(rownames(Mut)[y1])
    B=Mut.sam
    
    #keep the interactions with at least one gene expressed differently
    U.edge=Dysnet
    Mutdys=c()
    for(i in 1:dim(U.edge)[1]){
      xa=which(sigwilcox[,"signature"]==U.edge$GeneA[i])
      xb=which(sigwilcox[,"signature"]==U.edge$GeneB[i])
      if(length(xa)>0|length(xb)>0){
        Mutdys=rbind(Mutdys,U.edge[i,])
      }
    }
    
    if(is.null(Mutdys)){ 
      #
      return(NULL)
    }
    
    #evaluation of co-occurrence
    MutEdgeP <- tapply(Mutdys$Sam.ID,paste(Mutdys$GeneA,Mutdys$GeneB,sep=";"),function(samples){
      Dys.sam=unique(samples)
      A=Dys.sam
      #B=Mut.sam
      #Simulation p-value
      sim=unlist(lapply(1:n.sim,
                        function(i){AA=sample(1:N11,length(A));BB=sample(1:N22,length(B));return(sum(AA %in% BB))}))
      P=length(which(sim>=length(intersect(A,B))))/n.sim
      
      WKK <- data.frame(mutgene,length(B),length(A),length(intersect(A,B)),P)
      colnames(WKK) <- c("MutGene","MutSamples","DysSamples","CooccuSamples","P")
      return(WKK)
      
    })
    
    FF.net0 <- do.call(rbind,MutEdgeP)
    FF.net1 <- cbind(rownames(FF.net0),FF.net0)
    colnames(FF.net1)[1] <- c("GeneA-B")
    
    return(FF.net1)
  })
  
  FF.net <-do.call(rbind,FF.nets)
  
  Sig=which(FF.net$P<alpha)
  MutPPI=FF.net[Sig,]
  
  res <- list(MutPPI=MutPPI,FF.net=FF.net)
  return(res)
}

