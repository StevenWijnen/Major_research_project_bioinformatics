ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"
library(BSgenome)
library(dplyr)
library(org.Hs.eg.db)
library(tidyr)
library(VariantAnnotation)
library(ref_genome, character.only = TRUE)
library(purrr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


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

    patient_data.2 <- patient_data %>% merge(extensive_patient_data, by="Biomaterial_Id")
    
    #Add group to the dataframe
    patient_data.2$Group <- ifelse(patient_data.2$short %in% c('AML', "bALL"), "low_IRR", "high_IRR")
}

