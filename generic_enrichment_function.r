# Enrichment function
GSHG_custom <- function(allModules, # modules to be analyzed
                 mySets, # a list of all genesets
                 allgenes # all genes in dataset being analyzed
                 ){
  fisherTest <- function(set,mod,all){
    total.shared <- length(intersect(all,set))
    shared.in.mod <- length(intersect(mod,set))
    shared.out.mod <- total.shared-shared.in.mod
    in.mod.not.shared <- length(mod)-shared.in.mod
    out.mod.not.shared <- length(all)-length(mod)-shared.out.mod
    fisher.test(matrix(c(shared.in.mod,in.mod.not.shared,shared.out.mod,out.mod.not.shared),ncol=2),
                alternative="greater")$p.val
  }  
  
  allnetGenes <- unique(unlist(allModules))
  checkOverlap=function(x,y){
    length(intersect(x,y))
  }
  
  intsct <- sapply(mySets,checkOverlap,allnetGenes)
  mySets <- mySets[intsct > 0]
 # cat("Kept",length(mySets), "sets")
  
  mySetNames <- names(mySets)
  
  GSHGresults=matrix(nrow=length(mySets),
                     ncol=length(allModules),
                     data=-8888)
  
  for(i in c(1:length(allModules))){
    allModGenes <- unlist(allModules[[i]])
    GSHGresults[,i] <- unlist(lapply(mySets,
                                     function(aSet)
                                       fisherTest(aSet,allModGenes,all=allgenes)))
  #  print(paste("Finished module",i))
  }

  for(pla in seq_along(mySets)){
    fisherTest(mySets[[pla]], allModGenes, all=allgenes)
  }

  colnames(GSHGresults) <- names(allModules)
  datout <- data.frame("sets"=mySetNames,
                          GSHGresults)
  return(datout)
}