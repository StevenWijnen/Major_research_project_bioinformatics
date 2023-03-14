#This script also performs a label permutation but here the expensive kfold calculation is ommited 
#Since we are not interessted in kappa scores but the calculated cat scores for each random model

library(clusterProfiler)
library(dplyr)
library(pryr)
library(foreach)
library(doParallel)
library(doRNG)
library(sda)

wd <- getwd()
source(paste(wd, "RNA/Scripts/LoadData.R", sep="/"))
#Get the kegg pathways for hsa

kegg <- download_KEGG("hsa")

#INPUT: data -> dataframe: The dataframe for which you want to perform permutation test
#       group-> string: either cat.os+tall for TARGET or cat.os+tall for PMC; This represents the group within the dataframe that needs to be used to calculate the T^2 
#       dataset -> string: helper to write file 
#       output_dir -> string: path to output directory
main <- function(data, group, dataset="", output_dir="") {
  X <- data[["X"]] 

  #We needed a binary phenotype so we are using Group models
  y <- data[["Y"]]$Group

  ans <- parallel_func(X,y, group)

  
  write.csv(ans, file= paste(
    paste(wd, output_dir, sep="/") , 
    paste(dataset, "_kegg_permutation_20000", sep=""), 
    sep="/" ))  
  
  return(ans) 
}

#Again we use multi threading
parallel_func <- function(X, y, group) {
  print(paste("Memory in use", mem_used()))
  start_time <- Sys.time()
  
  #Create a cluster
  cl <- makeCluster(47, methods = FALSE, type="FORK")
  
  #export functions nessecarry in the single threads
  #For PMC HPC 47 cores were used
  clusterExport(detectCores()-1, c("tsquare_kegg"))
  registerDoParallel(cl)
  registerDoRNG(seed = 1616211)
  
  #There are 352 kegg pathways for which we need to report the T^2 value of the random model
  #Create empty matrix to speed up parallelization
  temp_df <- matrix(data = NA, nrow = 352, ncol = 20003)
  temp_df <- foreach(i=1:20000, .combine=cbind) %dopar% {
      print(paste("number:", i))
      yhat <- sample(y)
    
    ra <- sda.ranking(X, yhat, fdr=TRUE, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy")
    nFeatures <- sum(ra[,"lfdr"] <= 1)
    
    selVars = ra[,"idx"][1:nFeatures]
    
    df <- ra[,] %>% as.data.frame() %>% tibble::rownames_to_column(var="Gene")

    real <- tsquare_kegg(df, group)
    #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    real 
  }

  stopCluster(cl)
  end_time <- Sys.time()
  
 print(paste("time spent...:", (end_time - start_time)))
  
  return(temp_df)
}

#Calculate the T^2 value of a pathway with the random permuted labels
tsquare_kegg <- function(df, group, dir="") {
  library(org.Hs.eg.db)
   
  #some genes are then duplicated (only 1 or 2) so only keep the first one
  #apperently -c(0) returns empty dataframe, so first check if there are duplcites
  if(length(which(duplicated(df$Gene))) > 0) {   
     df <- df[-which(duplicated(df$Gene)),]
  } 
  #Map ensembl ids to entrezids used in KEGG
  mapper <- bitr(df$Gene, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
 
  list_T2<- c()
  #for every kegg pathway...
  for(pathway in unique(kegg$KEGGPATHID2EXTID$from)) {
    #Get the entrez ids of the pathway
    entrez_ids_of_pathway <- kegg$KEGGPATHID2EXTID[kegg$KEGGPATHID2EXTI$from == pathway,"to"]
    
    #map entrez ids to ensembl ids
    ensembl_ids <- mapper[mapper$ENTREZID %in% entrez_ids_of_pathway, "ENSEMBL"]
    ensembl_ids <- ensembl_ids[!duplicated(ensembl_ids)]
    
    #Calculate T^2
    T2 <- sum(df[df$Gene %in% ensembl_ids, group]^2 )
    gene_specific_score <- data.frame("gene" = ensembl_ids, 
                                      "value" = sign(df[df$Gene %in% ensembl_ids, group]) * sqrt(T2)
    )
    
    directed_scores <- c()
    #Add regulation directions if supplied
    if(dir != "") {
      genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group, dir)]
      directed_scores <- abs(genes_of_intrest[,group]) * genes_of_intrest[, dir]
      names(directed_scores) <- genes_of_intrest$Gene
    }
    
    #Save genes found and T^2 value 
    genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group)]
    list_T2[pathway] <-   T2
  }

   return(list_T2)
}



target_data <- load_TARGET_data()
main(target_data, "cat.os+tall", "TARGET", "RNA/Permutation_Tests/output")


#TODO uncomment for PMC analysis, data can be requested
#pmc_data <- load_PMC_data()
#colnames(pmc_data[["X"]]) <- substr(colnames(pmc_data[["X"]]),1,15)
#main(pmc_data, "cat.os+tall", "PMC", "RNA/Permutation_Tests/output")

