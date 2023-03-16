#Script for the filtering of the vcf files

library(BSgenome)
library(dplyr)
library(org.Hs.eg.db)
library(tidyr)
library(VariantAnnotation)
library(purrr)

#PMC data works with hg38 reference genome 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


#VCF file tumor and normal sample need a dp of >=10
#Input for these filters is always RowRange object from the vcf file
dp_filter <- function(x) {
  dp <- geno(x)$DP
  normalDP <- dp[,1]
  tumorDP <- dp[,2]
  
  b <- (normalDP >= 10) & (tumorDP >= 10) 

  return(b)
}

#The code for the vaf filter needed some additional steps, this was needed since the order of normal tumor sample in the vcf files was not 
#consistent in the files. That's why when looping over the tumor_ids the vaf filter finciton is created with the correseponding tumor id in place
#This id is then used to check if the vaf of the tumor sample is high enough (>=0.1)

wrapperFunction <- function(tumor_id) {
  vaf_filter <-  function(x) {
    af <- geno(x)$AF
    af <- af[,tumor_id]
    b <- sapply(af, function(l) {
      if(length(l) == 1 &&l >= 0.1) {
        return(TRUE)
      }
      return(FALSE)
    })

    return(b)
  }
  return(vaf_filter)
}

#Quality filter MMQ filter for both tumor and normal sample; here the order of tumor normal doesn't matter since we are checking both 
#And both need to be >= 60
mmq_filter <- function(x) {
  mmq <- info(x)$MMQ
  tumor_mmq<-sapply(mmq, function(ls) {ls[2]})
  normal_mmq<-sapply(mmq, function(ls) {ls[1]})

  b1 <- tumor_mmq >= 60
  b2 <- normal_mmq >= 60
  b <- b1 | b2
  
  return(b)
}

#Loop over the directory containing vcf and tbi files and filter the vcf file
#Than write the filtered file to the right location
filter_vcf_files <- function(output_dir="") {
    #Get the files that needs to be filtered
    patient_data <- read.csv("WGS/data/data_request_01_patient_data.csv",  sep=";")
    
    #Parse all vcf files for which we have at least tumor WGS data;
    #Obtain there ids
    tumor_ids <- patient_data %>% 
        filter(X01_Biomaterial_type == "DNA") %>%   
        filter(X02_Disease_status == "tumor") %>% 
        dplyr::select(Biomaterial_Id) %>% 
        unlist()

    #Here we only get the file names
    files <- list.files("WGS/data/unfiltered_vcf/", pattern="*.vcf.gz$")
    #Here we append the path to the file name
    files_locations <- paste("WGS/data/unfiltered_vcf/", files, sep="/")

    for(id in tumor_ids) {
        if(any(grepl(id, files_locations))) {
            index <- which(grepl(id, files_locations))
            #Create the vaf filter with correct tumor id
            vaf_filter <- wrapperFunction(id)
            
            #Set all the filters
            filters <- FilterRules(list(DP=dp_filter, VAF=vaf_filter, MMQ= mmq_filter))
            destination.file <- paste(output_dir, substr(files[index],1, nchar(files[index])-3 ), sep="/")
            
            #Perform the actual filtering
            filterVcf(files_locations[index], "hg38", destination.file, filters=filters, verbose=TRUE)
        }
    }
}   

filter_vcf_files("WGS/output_filtered_vcf")





