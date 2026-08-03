library(data.table)
library(tidyverse)
library(qs)

## Some exploratory analyses calculating similarities between modules
# Calculate jaccard similarities between two networks
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

path1 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/Bicor-no_TO_signum0.287_minSize5_merge_ME_0.85_17193/kME_table_.csv"))
path1m <- tapply(path1[,1], path1[,2], list) # 48 mods
path1m <- path1m[-which(names(path1m)=="turquoise")]
path2 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/Bicor-no_TO_signum0.38_minSize3_merge_ME_0.85_17193/kME_table_.csv"))
path2m <- tapply(path2[,1], path2[,2], list) # 18 mods
path2m <- path2m[-which(names(path2m)=="turquoise")]

jsim <- lapply(path1m, function(x) lapply(path2m, function(y) jaccard(x,y)) %>% unlist)
lapply(jsim,max) %>% unlist %>% sort(decreasing=T) # 10 mods out of 47 greater than 0 -> 10 of the bulkMinusLeinpb mods 
lapply(jsim, function(x){
  if(mean(x)==0){
    return("none")
  } else {
    return(names(x)[which.max(x)])
  }
}) %>% unlist %>% unique # brown blue yellow <- these mods are from the megaMinusSCVI network that 

jsim2 <- lapply(path2m, function(x) lapply(path1m, function(y) jaccard(x,y)) %>% unlist)
lapply(jsim2,max) %>% unlist %>% sort(decreasing=T) # 15 mods out of 17 greater than 0
lapply(jsim2, function(x){
  if(mean(x)==0){
    return("none")
  } else {
    return(names(x)[which.max(x)])
  }
}) %>% unlist %>% unique # blue brown red


