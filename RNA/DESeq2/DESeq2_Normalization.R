#NOTE, data for this script is only available upon request and not present in this repository

library(DESeq2)
library(dplyr)
library(sva)
library(readr)

#Not uploaded since it contains sensitive information
#Can be retrieved upon request
source( paste(wd, "RNA/DESeq2/weird_samples.R", sep="/"))


main <- function(output_dir="") {
  #Read patient data from PMC
  df_patient <- read_patient_data("/hpc/pmc_vanboxtel/personal/swijnen/SDA/data/data/PMC/RNA_patients.csv")
  #Read raw counts
  df_cts_PMC <- read_count_files(dir="/hpc/pmc_vanboxtel/projects/SexDifferences/1_Input/RNAseq_PMC")
  #Perfrom DESEq2 analysis
  desq_data <- desq_pipe(df_cts_PMC, df_patient)

  #Save the patient data, normalised count and vsd matrix
  saveRDS(df_patient, file= paste(output_dir, "df_patients.RData", sep="/"))
  saveRDS(desq_data[[1]], file = paste(output_dir, "normalised_counts.RData", sep="/"))
  saveRDS(desq_data[[2]], file = paste(output_dir, "vsd.RDta", sep="/"))
  
  print("finished")
}

#Input: raw counts matrix; patient data
#Output: list containing vst normailized counts and DESeq2 normalized counts 
desq_pipe <- function(df_cts_PMC, PMC_patients) {
  cts <- df_cts_PMC
  rownames(cts) <- rownames(df_cts_PMC)
  coldata <- PMC_patients[PMC_patients$Biomaterial_Id %in% colnames(cts),c("Biomaterial_Id", "Gender", "Tumor_type", "Tissue", "Disease")]  %>% 
    arrange(Biomaterial_Id)
    
  #only use main groups
  coldata <- coldata %>% dplyr::filter(Disease %in% c("bALL", "tALL", "AML","Burkitt"))
  
  #Biomaterial ids from older pipeline version, containing a batch effect which needs to be removed.
  #These different cancer types cluster together in the PCA if not corrected for..
  weird_ones <- get_weird_samples()
  coldata <- coldata %>% mutate(batch = ifelse(Biomaterial_Id %in% weird_ones, "b2", "b1"))
  cts<- cts[, coldata$Biomaterial_Id]

  #Make Gender column a categroical variable   
  coldata$Gender <- as.factor(coldata$Gender)
  rownames(coldata) <- coldata$Biomaterial_Id
  
  #If the rownames (patient_ids) of the count matrix and column matrix are not the same throw an error.
  if(!(all(rownames(coldata) %in% colnames(cts)) & all(rownames(coldata) == colnames(cts)))) {
    stop("ERROR: Count matrix data is not the same as column data")
  } else {
    #Remove btach effect with combat 
    #Create design matrix for covariates we want to keep; here gender
    design0 <- model.matrix(~ coldata$Gender)
  #If biomaterial id belongs to the old batch set batch to b2
	coldata<-coldata %>% mutate(batch = ifelse(Biomaterial_Id %in% weird_ones, "b2", "b1"))
	#Use batch correction tool on the original cts matrix. Here we design the matrix such that  Gender and Disease such that effects aren't omitted by the tool
  corrected_cts <- ComBat_seq(data.matrix(cts), coldata$batch, group=coldata$Disease, covar_mod=design0)   
    
   #Perform DESeq2 analysis  
	dds <- DESeqDataSetFromMatrix(countData = corrected_cts,
                                  colData = coldata,
                                  design = ~Disease)
        #Calculate size factors
   dds <- estimateSizeFactors(dds)
    sizeFactors(dds)
    #counts scales the counts to library size to account for differences in sequencing depth across samples.
    normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
    rownames(normalized_counts) <- rownames(df_cts_PMC)
    
    #Normalization used for PCA analysis
    rld <- vst(dds, blind=TRUE)
    
    return(list(normalized_counts, rld))
    
  }
  return(0)
}

#Load count data
#Input: Directory contaitng all count files you want to analyse
#Output: matrix containg counts per gene per sample (patient)
read_count_files <- function(dir="") {
  #Extract file names from directory
  count_files_names <- list.files(dir, pattern="*.txt")
  #Read first file name to create initial dataframe
  df_cts_PMC <- read_delim(paste(dir,count_files_names[1], sep="/"), comment = "#", col_names=c("ID",	"Counts",	"CPM"	,"FPKM",	"Chr",	"Start",	"End"	,"Strand"	,"Length"	,"GeneID"	,"GeneName",	"TranscriptID"), delim="\t")
  #Only keep gene id and raw counts
  df_cts_PMC <-  df_cts_PMC[,c("ID", "Counts")]
  
  colnames(df_cts_PMC)[2] <- sub("_.*","", count_files_names[1])
  
  #Loop over all files and append to the count matrix
  for(file in count_files_names[2:length(count_files_names)]) {
    #Extract biomaterial Id from the file name
    bio_mat_id <- sub("_.*","", file)
    #I've tested the double once, either one of the files was smaller in lines, hence the next argument 
    #otherwise the two files were interchangable in the PCA plot, so just keep the first 
    if(!(bio_mat_id %in% colnames(df_cts_PMC))) {
      df <- read_delim(paste(dir,file, sep="/"), comment="#", col_names=c("ID",	"Counts",	"CPM"	,"FPKM",	"Chr",	"Start",	"End"	,"Strand"	,"Length"	,"GeneID"	,"GeneName",	"TranscriptID"), delim="\t")
      if(dim(df)[1] < 60000) {
        #Some files were double of which one contained fewer genes, skip these files
        print(paste("not that much genes ", file, as.character(dim(df)[1])))
        next
      }
      
      df <- df[,c("ID", "Counts")]
      colnames(df)[2] <- bio_mat_id
      #Append to matrix
      df_cts_PMC <- merge(df_cts_PMC,df, by="ID")
    } else {
      print(paste("double" ,file))
    }
  }
  #Set rownames for further analysis
  rownames(df_cts_PMC) <- df_cts_PMC$ID
  df_cts_PMC <- df_cts_PMC[,-1]
  
  return(df_cts_PMC)
}

#Input: file location of patient data
#Outpur: dataframe containing patient information
read_patient_data <- function(file="") {
  PMC_patients <- read_csv(file)[,-1]
  PMC_patients$Tumor_type <- as.factor(PMC_patients$Tumor_type)
  #The patient file contains long names for disease, map these to shorter names for convenience 
  mapper <- data.frame("Tumor_type"= c(
  "Acute leukemia, NOS"                 ,                                           
  "Acute myeloid leukaemia, t(8;21)(q22;q22)"       ,                               
  "Acute myeloid leukemia with abnormal marrow eosinophils (includes all variants)",
  "Acute myeloid leukemia, minimal differentiation"                        ,        
  "Acute myeloid leukemia, NOS (FAB or WHO type not specified, see also M-9930/3)" ,
  "Acute myelomonocytic leukemia"                   ,                               
  "Acute promyelocytic leukaemia, t(15;17)(q22;q11-12)"       ,                     
  "Aggressive NK-cell leukaemia"       ,                                            
  "B lymphoblastic leukemia/lymphoma, NOS"           ,                              
  "Burkitt cell leukemia (see also M-9687/3)"              ,                        
  "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)"       ,        
  "Juvenile myelomonocytic leukemia"  ,                                             
  "Mixed phenotype acute leukemia, B/myeloid, NOS"     ,                            
  "Myeloid leukemia, NOS"  ,                                                        
  "Precursor B-cell lymphoblastic leukemia (see also M-9728/3)"     ,               
  "Precursor T-cell lymphoblastic leukemia (see also M-9729/3)" 
  ), 
    
  Disease=c(
    rep("AML",5),
    "Myelomonocytic",
    "promyelocytic",
    "NKL",
    "blymphoma", 
    "Burkitt_leukemia",
    "Burkitt",
    "myelomonocytic",
    "Mixed",
    "AML", 
    "bALL", 
    "tALL"
  ))

  #Sort by Biomaterial_Id such that the order of patient data is inline with the rownames of the cts matrix
  PMC_patients <- merge(PMC_patients, mapper, by="Tumor_type") %>% arrange(Biomaterial_Id)
  
  return(PMC_patients) 
}

main()