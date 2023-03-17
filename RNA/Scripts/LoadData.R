#Script to load the normalized count data for the TARGET and the PMC dataset

library(readr)
library(dplyr)
library(tidyr)

wd<- getwd()
#setwd(paste(wd, "RNA/", sep="/"))
source(paste(wd, "RNA/Scripts/blacklist_pmc_ids.R", sep="/"))

load_TARGET_data <- function() {
  dir <- paste(wd, "RNA/data/TARGET/", sep="/")
  #Read patient files
  tALL_patients <- (read_csv(paste(dir, "tALL_patients.csv", sep=""), show_col_types=FALSE))[,3:5]
  bALL_patients <- (read_csv(paste(dir,"bALL_patients.csv", sep=""), show_col_types=FALSE))[,3:5]
  AML_patients <- (read_csv(paste(dir,"AML_patients.csv",sep=""), show_col_types=FALSE))[,3:5]
  OS_patients <- (read_csv(paste(dir,"OS_patients.csv", sep=""), show_col_types=FALSE))[,3:5]
  NBL_patients <- (read_csv(paste(dir,"NBL_patients.csv", sep=""), show_col_types=FALSE))[,3:5]
  
  #Read TPM normalized count files
  tALL_tpm <- (read_csv(paste(dir,"tALL_tpm.csv", sep=""), show_col_types=FALSE))[,-1]
  bALL_tpm <- (read_csv(paste(dir,"bALL_tpm.csv", sep=""), show_col_types=FALSE))[,-1]
  AML_tpm <- (read_csv(paste(dir,"AML_tpm.csv", sep=""), show_col_types=FALSE))[,-1]
  OS_tpm <- (read_csv(paste(dir,"OS_tpm.csv", sep=""), show_col_types=FALSE))[,-1]
  NBL_tpm <- (read_csv(paste(dir,"NBL_tpm.csv", sep=""), show_col_types=FALSE))[,-1]
  
  #Gene mapper, for gene names to ENSEMBL ids as provided by TARGET 
  TARGET_mapper <- read_csv(paste(dir,"TARGET_gene_mapper.csv", sep=""))[,c(2,3)]
  #Only tALL file was in gene names
  tALL_tpm<- tALL_tpm %>% merge(TARGET_mapper, by.x="Gene", by.y="hgnc_symbol") %>% dplyr::select(-Gene)
  colnames(tALL_tpm)[dim(tALL_tpm)[2]] <- "Gene"
  
  #Merge the different cancer type count files, and put in wide format.
  combined_tpm <- tALL_tpm %>%
    merge(bALL_tpm, by="Gene") %>%
    merge(AML_tpm, by="Gene") %>%
    merge(OS_tpm, by="Gene") %>%
    merge(NBL_tpm, by="Gene") %>%
    pivot_longer(cols=-Gene, names_to="patient_id",values_to='tpm') %>%
    pivot_wider(names_from= Gene, values_from = "tpm")
  
  #Combine all patient data and provide grouping 
  Y <- rbind(
    tALL_patients, 
    bALL_patients,
    AML_patients, 
    OS_patients, 
    NBL_patients) %>% 
    mutate(Group = ifelse(Disease == "AML" | Disease == "NBL" | Disease=="bALL", "aml+nbl+ball", "os+tall") ) 
  
  #Sanity cehck if all patient ids are in both dataframes the same
  if(all(combined_tpm$patient_id == Y$patient_id)) {
    #Add the group_gender column by pating group and sex
    Y <- Y %>% mutate(Group_Gender = paste(Group, Gender, sep="_"))
    X <- make_matrix(combined_tpm)
    return(list("Y"=Y,"X"= as.matrix(X)))
  } else {
    print("Not all ids are in cts and patient file; TARGET")
    stop()
  }
}

#PMC data is confidential and can be obtained upon request
load_PMC_data <- function() {
  #use combat batch corrected normalized counts as calculated by DESeq2_Normalization script.
  dir <- paste(wd, "RNA/data/PMC/", sep="/")

  df_patients <-readRDS(paste(dir, "df_patients.RData", sep=""))
  normalized_counts <- readRDS(paste(dir,"normalised_counts.RData", sep=""))
  
   # filter out one blacklisted PMC patient 
  df_patients <- remove_id(df_patients)
  df_cts <- normalized_counts %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var="GeneId") %>% 
    pivot_longer(cols=-GeneId , names_to="patient_id", values_to="normalised_count") %>% 
    pivot_wider(names_from = "GeneId", values_from = "normalised_count")
  


  #test commit 
  coldata_total <- df_patients[df_patients$Biomaterial_Id %in% df_cts$patient_id,c("Biomaterial_Id", "Gender", "Tumor_type", "Tissue", "Disease")]  %>% 
    arrange(Biomaterial_Id)
  
  #Some tissue types are missing, provide string
  coldata_total$Tissue[is.na(coldata_total$Tissue)] <- "Tissue missing"
  #We are only intreseted in bALL, tALL, AML and Burkitt samples
  coldata_total_main_groups <- coldata_total %>% dplyr::filter(Disease %in% c("bALL", "tALL", "AML","Burkitt"))
  #Filter out patients wich are not from the cancer types which we are analysing
  df_cts <- df_cts %>% dplyr::filter(patient_id %in% coldata_total_main_groups$Biomaterial_Id)
  
  #Sanity check if the dataframes are in the same order
  if(all((coldata_total_main_groups$Biomaterial_Id == df_cts$patient_id))) {
    #Create grouping
    Y<- coldata_total_main_groups %>% 
      mutate(Group = ifelse(coldata_total_main_groups$Disease=="AML" | coldata_total_main_groups$Disease =="bALL", "aml+ball", "tall+b"))
    #Create group_gender column
    Y <- Y %>% 
      mutate(Group_Gender = paste(Group, Gender, sep="_"))
    
    X <- make_matrix(df_cts)    
    
    return(list("Y"=Y, "X"=X))
  }
  else {
    print("patient data and cts are different identifiers; PMC")
    stop()
  }
}

make_matrix <- function(df) {
  #Remove zero variance varibles since we cannot learn anything from them.
  #Removing them will make our consecutive analysis faster.
  X <- df[,-1]
  X <- as.matrix(X[, -caret::nearZeroVar(X)])
  rownames(X) <- df$patient_id
  
  
  return(X)
}


#load_PMC_data()
#load_TARGET_data()
