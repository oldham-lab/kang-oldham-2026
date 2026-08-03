# GSEA related functions

library(GSEABase)
library(limma)
library(qvalue)
library(ellipse)
library(WGCNA)
library(flashClust)
library(data.table)

GSHG=function(allModules,
              allgenes,
              set_list, # a list of genesets with names 
              file_desc,
              save_dir,
              save_name
){
  

  mySetNames <- names(set_list)
  mySets <- set_list
  
  mySets=lapply(mySets,unlist)
  mySets=lapply(mySets,as.character)
  mySets=lapply(mySets,toupper)
  
  AllModGenes=unique(unlist(allModules))
  
  checkOverlap=function(x,y){
    length(intersect(x,y))
  }
  
  intsct=sapply(mySets,checkOverlap,AllModGenes)
  keptSets=intsct>0
  mySets=mySets[keptSets]
  mySetNames=mySetNames[keptSets]
  file_desc = file_desc[keptSets]
  
  GSHGresults=matrix(nrow=length(mySets),ncol=length(allModules),data=-8888)
  
  allgenes=allgenes[!is.na(allgenes)&allgenes!=""]
  allgenes=gsub(" ","",allgenes)
  
  fisherTest=function(set,mod,all){
    total.shared=length(intersect(all,set))
    shared.in.mod=length(intersect(mod,set))
    shared.out.mod=total.shared-shared.in.mod
    in.mod.not.shared=length(mod)-shared.in.mod
    out.mod.not.shared=length(all)-length(mod)-shared.out.mod
    fisher.test(matrix(c(shared.in.mod,in.mod.not.shared,shared.out.mod,out.mod.not.shared),ncol=2),alternative="greater")$p.val
  }
  
  for(i in c(1:length(allModules))){
    allModGenes=unlist(allModules[[i]])
    GSHGresults[,i]=unlist(lapply(mySets,function(aSet){fisherTest(aSet,allModGenes,all=allgenes)}))
   # print(paste("Finished module",i))
  }
  
  GSHGresults[abs(GSHGresults)==0]=1e-300
  colnames(GSHGresults)=names(allModules)
  
  datout=data.frame("SetID" = mySetNames,
                    "SetName" = file_desc,
                    "SetSize" = unlist(lapply(mySets, length)),
                    as.data.frame(GSHGresults))
  #fwrite(datout, file = paste0(save_dir, "/", save_name, ".csv"))
  return(datout)
} ## End of function

BroadGSHG <- function(allModules,
                      allgenes,
                      allSets,
                      save_dir,
                      save_name){

 
	AllModGenes=unique(unlist(allModules))
	
	checkOverlap=function(x,y){
		length(intersect(x,y))
	}
	
	intsct=sapply(geneIds(allSets),checkOverlap,AllModGenes)
	keptSets=intsct>0
	allSets=allSets[keptSets]
	setIDs=sapply(allSets,setIdentifier)
	collType=lapply(allSets, collectionType)
	catType=sapply(collType, bcCategory)
	setGenes=sapply(allSets,geneIds)
	setSize=sapply(setGenes,length)
	species=sapply(allSets,organism)
	descriptions=sapply(allSets,description)
	pubmed=as.character(sapply(allSets,pubMedIds))
	pubmed[pubmed==""]=NA
	 
	GSHGresults=matrix(nrow=length(allSets),ncol=length(allModules),data=-8888)
	
	allgenes=allgenes[!is.na(allgenes)&allgenes!=""]
	allgenes=gsub(" ","",allgenes)
	
	fisherTest=function(set,mod,all){
		
		totalshared=length(intersect(all,set))
		modshared=length(intersect(mod,set))
		nonmodshared=totalshared-modshared
		modnonshared=length(mod)-modshared
		nonmodnonshared=length(all)-length(mod)-nonmodshared
		fisher.test(matrix(c(modshared,modnonshared,nonmodshared,nonmodnonshared),ncol=2),alternative="greater")$p.val
		
	}
	
	for(i in c(1:length(allModules))){
		allModGenes=unlist(allModules[[i]])
		GSHGresults[,i]=unlist(parallel::mclapply(allSets,function(aSet){fisherTest(geneIds(aSet),allModGenes,all=allgenes)}, mc.cores=8))
		print(paste("Finished module",i))
	}
	
	GSHGresults[abs(GSHGresults)==0]=1e-300
	colnames(GSHGresults)=names(allModules)
  datout=data.frame("SetID" = setIDs,
                    "SetName" = names(allSets),
                    "SetSize" = setSize,
                    as.data.frame(GSHGresults))
 # fwrite(datout, file = paste0(save_dir, "/", save_name, ".csv"))
  return(datout)
}

run_gsea_for_proj <- function(allModules,
                              set_list=NULL,
                              file_desc=NULL,
                              broad = T,
                              save_dir
                              ){
  allgenes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[,2]

  if(is.null(set_list)){
    
    # load Oldham genesets
    set_dir <- system.file("extdata", "GeneSets_consensus", package = "CoPA")
    gsea_legend <- read.csv(paste0(set_dir,"/ConsensusGeneSets.csv"))
    files <- list.files(set_dir, full.names=T)[grep("MOSET", list.files(set_dir))]
    file_names <- list.files(set_dir)[grep("MOSET", list.files(set_dir))]
    set_list <- list()
    file_desc <- c()
    for(l in 1:length(files)){
      set_list[[l]] <- read.csv(files[l], header=F)
      names(set_list)[l] <- gsub(".csv","",file_names[l])
      file_desc[l] <- as.character(gsea_legend[gsea_legend[,1] == names(set_list)[[l]],2])
    }
  }

  if(broad){
    broadSets <- getBroadSets(Sys.getenv("MSIGDB_XML", "/home/gugene/code/git/COPA/data/msigdb_v7.4.xml"))
    cat("Running GSEA using Broad genesets...\n")
    BroadGSHG(allModules = allModules,
              allgenes = allgenes,
              allSets = broadSets,
              save_dir = save_dir,
              save_name = "gsea_Broad")
  } else {
    cat("Running GSEA using user input genesets...\n")
    GSHG(allModules = allModules,
         allgenes = allgenes,
         set_list = set_list,
         file_desc = file_desc,
         save_dir = save_dir,
         save_name = "gsea_userInput")
  }
  
}