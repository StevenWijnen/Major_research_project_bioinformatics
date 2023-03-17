library(dplyr)
library(org.Hs.eg.db)
library(tidyr)
library(VariantAnnotation)
library(purrr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(ggpubr)

setwd("/Users/stevenwijnen/surfdrive/PMC/Major_research_project_bioinformatics/")
#This file contains sensitive PMC ids, only upon requst can these ids be shared
source("WGS/wrong_samples.R")

#List of all found genes with their correseponding KEGG pathway
load("WGS/data/kegg_t2.RData")

#Load the reference genome of known genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Code to load in the patient data
load_patient_data <- function() {
  #Load in the filtered patients 
  patient_data <- read_csv("WGS/data/WGS_patients.csv")
  
  #Code to map large disease names to smaller once for conveinince 
  tumor_types <- c(
    "Acute myeloid leukemia, NOS (FAB or WHO type not specified, see also M-9930/3)",
    "Acute myeloid leukaemia, t(8;21)(q22;q22)",
    "Precursor B-cell lymphoblastic leukemia (see also M-9728/3)",
    "Precursor T-cell lymphoblastic leukemia (see also M-9729/3)",
    "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)"
  )
  #Filter only the cancer types of intrest 
  patient_data <- patient_data %>%filter(Tumor_type %in% tumor_types)
  mapper <- data.frame("Tumor_type"=levels(as.factor((patient_data$Tumor_type))), short=c(
    rep("AML",2),
    "BL",
    "bALL",
    "tALL"
  ))
  
  patient_data <- merge(patient_data, mapper, by="Tumor_type")
  
  patient_data<- patient_data %>% filter(Disease_status == "tumor")
  
  #Some meta daa was lost in the filtered file, add the meta data to the filtered patients
  extensive_patient_data <- read.csv("WGS/data/data_request_01_patient_data.csv",  sep=";")
  
  patient_data_extensive <- patient_data %>% merge(extensive_patient_data, by="Biomaterial_Id")
  
  #Add group to the dataframe
  patient_data_extensive$Group <- ifelse(patient_data_extensive$short %in% c('AML', "bALL"), "low_IRR", "high_IRR")
  
  return(patient_data_extensive)
}

#Code that checks all the coding mutations in the samples
overall_coding_mutations <- function() {
  #Get the filtered vcf files
  files <- list.files("WGS/output_filtered_vcf", pattern=".vcf")
  files <- paste("WGS/output_filtered_vcf",files,sep="/")

    #Create empty dataframe
  df_coding_variant_mutations <- data.frame()
  
  patient_data <- load_patient_data()
  
  for(file in files) {
    #Get the id of the sample
    id <- substr(file, 34,44)
    
    if(id %in% patient_data$Biomaterial_Id) {
      print(id)
      vcf <- readVcf(file, "hg38")
      available_seqlevels<-seqlevels(vcf)[which(seqlevels(vcf) %in% seqlevels(txdb))]
      vcf<-keepSeqlevels(vcf, available_seqlevels)
      rd <- rowRanges(vcf)
      genes_queried_intresting <- locateVariants(rd, txdb, CodingVariants())
      
      Gender <- patient_data[patient_data$Biomaterial_Id == id, "X03_Sex"]
      Group <- patient_data[patient_data$Biomaterial_Id == id, "Group"]
      Disease <- patient_data[patient_data$Biomaterial_Id == id, "X01_Tumor_type"]
      short <- patient_data[patient_data$Biomaterial_Id == id, "short"]
      df <- c("file" = file, "biomat_id"= id, "mutations"=length(genes_queried_intresting), "Gender"=Gender, "Group" = Group, "Group_gender" = paste(Group, Gender, sep="_"), "Disease"=Disease, "short"=short)
      df_coding_variant_mutations <- rbind(df_coding_variant_mutations, df)
    }
  }
  
  #Set the column names
  colnames(df_coding_variant_mutations) <- c("file" , "biomat_id",  "mutations", "Gender", "Group", "Group_gender", "Disease", "short")
  
  #Make the mutations numeric 
  df_coding_variant_mutations$mutations <- as.numeric(df_coding_variant_mutations$mutations)
  
  
  wrong_samples <- wrong_samples()
  
  df_coding_variant_mutations_filtered <-  df_coding_variant_mutations %>% filter(! biomat_id %in%  wrong_samples)
  
  return(df_coding_variant_mutations_filtered)
  
}

coding_mutations <- overall_coding_mutations()

ggplot(coding_mutations, aes(x = "", y =mutations, color=short, label=biomat_id))  +
  geom_boxplot(notch=T) +  
  xlab("Cancer type") + 
  ylab("Number of coding mutations") + 
  labs(color = "Cancer types") 


ggplot(coding_mutations, aes(x = "", y =mutations, color=Gender))  +
  geom_boxplot(notch=T) +  
  xlab("Gender") + 
  ylab("Number of coding mutations") + 
  scale_color_discrete(labels=c('Female', 'Male')) + 
  stat_compare_means(method = "wilcox.test")

summary(anova(lm(mutations ~ Gender + short, data=coding_mutations)))


genes_of_interest <- function() {
  #ENTREZ IDS of EXIT genes
  ATRX <- "546"
  CNKSR2 <- "22866"
  DDX3X <- "1654"
  KDM5C <- "8242"
  KDM6A <- "7403"
  MAGEC3 <- "139081"
  ZFX <- "7543"  
  USP9X <- "8239"
  
  EXIT_Genes <-  c(ATRX,CNKSR2, DDX3X, KDM5C, KDM6A, MAGEC3, ZFX, USP9X)
  
  intresting_genes <- list(gene_id=c(EXIT_Genes))
  
  KEGG_genes <- list()
  #Foreach pathway 
  for(name in names(all)) {
    #Get the genes corresponiding to the pathway 
    genes<- all[[name]]
    #Map ensembl ids to entrex ids
    mapper <- bitr(genes, "ENSEMBL", c("ENTREZID"), org.Hs.eg.db)
    
    KEGG_genes[[name]] <- mapper$ENTREZID
    intresting_genes[[1]] <- append(intresting_genes[[1]], KEGG_genes[[name]])
  }
  
  #Append the KEGG pathway genes to the XIT genes
  intresting_genes[[1]] <- intresting_genes[[1]][!duplicated(intresting_genes[[1]])]
  
  #Transform gene ids to GRANGES to perform the analysis
  granges_genes_of_interest <- mapIdsToRanges(txdb, keys =intresting_genes , type = "gene")
  
  return(list("granges" =  granges_genes_of_interest, "KEGG_genes" = KEGG_genes, "EXIT_genes" = EXIT_Genes, "intresting_genes"=intresting_genes))
}

map_mutations_to_genes <- function() {
  
  patient_data <- load_patient_data()
  ls <- genes_of_interest()
  granges <- ls[["granges"]]
  intresting_genes <- ls[["intresting_genes"]]
  
  #Get the filtered VCF files
  files <- list.files("WGS/output_filtered_vcf", pattern=".vcf")
  files <- paste("WGS/output_filtered_vcf",files,sep="/")
  
  #Create a dataframe with a patient id combined with each TXID of each gene of interest
  df_patients_with_intresting_gene_mutations_PMC <- as.data.frame(expand.grid("PMC.ID" = patient_data$Biomaterial_Id, "TXID" = intresting_genes[[1]]))
  #Set the mutation count for each patient-gene combination to 0 
  df_patients_with_intresting_gene_mutations_PMC$count <- 0
  
  used_ids <- c()
  #Get the genes of interest
  count <-0 
  for(file in files) {
    id <- substr(file, 34, 44)
    
    #id <- substr(file, 34,44)
    count <- count + 1
    print(count)
    if(id %in% patient_data$Biomaterial_Id) {
      print(id)
      used_ids <- c(used_ids, id)

      #Read the vcf file
      vcf <- readVcf(file, "hg38")
      
      #Get the available seqlevels for the taxonmy database
      available_seqlevels<-seqlevels(vcf)[which(seqlevels(vcf) %in% seqlevels(txdb))]
      #Keep only the seq levels which are avialble
      vcf<-keepSeqlevels(vcf, available_seqlevels)
      #Obtain the row ranges of the vcf 
      rd <- rowRanges(vcf)
      #Check if the patient has any coding mutations in one of the genes of interest
      genes_queried_intresting <- locateVariants(rd, granges, CodingVariants(), ignore.strand=TRUE)
      #If the patient has mutations...
      if(length(genes_queried_intresting) >0 ) {
        #Create a temporary datafframe
        df <- (data.frame(genes_queried_intresting))
        #Set the pmc id to the current tumor sample id
        df$PMC.ID<- id
        #Count the number of mutations per gene found
        ss <- df %>% group_by(PMC.ID, TXID) %>% summarise(i = n())
        
        #merge the dataframes by patient id and gene id now the new mutations of the current tumor id are appened to the dataframe
        #We will combine them later by adding the already found mutations to the new column
        #Therby adding only the new mutations of the patient while keeping the already found mutations
        df <- merge(df_patients_with_intresting_gene_mutations_PMC, ss, by=c("PMC.ID", "TXID"), all.x=TRUE)
        
        #Genes that were not found were set to NA remove these such that we can perform arthemetric operations again
        df[is.na(df$i),'i'] <- 0
        #Reduce the data frame back to three columns by adding the newly found mutations to the already known mutations
        
        df$count <- df$count + df$i
        
        #Remove the last column of the dataframe
        df_patients_with_intresting_gene_mutations_PMC<- df[,-4]
      }
    }
  }
  
  wrong_samples <- wrong_samples()
  df_patients_with_intresting_gene_mutations_PMC <- df_patients_with_intresting_gene_mutations_PMC %>% filter((!PMC.ID %in% wrong_samples), (PMC.ID %in% used_ids))
  

  return(df_patients_with_intresting_gene_mutations_PMC)
}


# #Know that we have the mutations perform the analysis
# mutations <- map_mutations_to_genes()
# interesting_genes <- genes_of_interest()

# KEGG_Genes <- interesting_genes[[2]]
# EXIT_Genes <- interesting_genes[[3]]

# patient_data <- load_patient_data()
# #EXIT gene analysis 

# #Get only the mutations in the EXIT genes
# df <- mutations %>%    
#   filter(TXID %in% EXIT_Genes) %>% 
#   group_by(PMC.ID) %>% 
#   summarise(total_mut = sum(count))  %>% 
#   distinct(PMC.ID, .keep_all=TRUE)  %>% 
#   merge(patient_data, by.x="PMC.ID", by.y="Biomaterial_Id") 

# #Create the boxplot of the EXIT genes
# ggplot(df, aes(x = "", y = total_mut, color=Gender))  +
#   geom_boxplot(notch=TRUE) + 
#   scale_color_discrete(labels = c("Female", "Male")) +
#   xlab("Gender") +
#   ylab("# Mutations EXIT genes")  + 
#   stat_compare_means(method = "wilcox.test")

# #Check if the significance still holds if we correct for the potential confounding effect disease
# m1 <- lm(total_mut ~ Gender + short , data=df)
# print(summary(m1))

# #For the KEGG pathway analysis only two pathways where significant plot these
# df1 <- mutations %>%
#   filter(TXID %in% KEGG_Genes[[9]]) %>%
#   group_by(PMC.ID) %>% summarise(total_mut = sum(count)) %>%
#   distinct(PMC.ID, .keep_all=TRUE)  %>%
#   merge(patient_data, by.x="PMC.ID", by.y="Biomaterial_Id")


# df2 <- mutations %>% 
#   filter(TXID %in% KEGG_Genes[[22]]) %>% 
#   group_by(PMC.ID) %>% summarise(total_mut = sum(count)) %>% 
#   distinct(PMC.ID, .keep_all=TRUE)  %>% 
#   merge(patient_data, by.x="PMC.ID", by.y="Biomaterial_Id")


# df1$Pathway = "PI3K-Akt signaling pathway"
# df2$Pathway = "Herpes simplex virus 1 infection"

# df_all <- rbind(df1,df2)

# ggplot(df_all, aes(x = Pathway, y = total_mut, color=Gender))  +
#   geom_boxplot(notch=TRUE) + stat_compare_means(method = "wilcox.test") +    
#   scale_color_discrete(labels = c("Female", "Male")) +
#   xlab("Pathway") +
#   ylab("# Mutations in pathway")
