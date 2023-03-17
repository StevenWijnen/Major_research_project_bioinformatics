library(readr)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
wd <- getwd()
source(paste(wd, "RNA/Scripts/LoadData.R", sep="/"))

levels_pmc <- c("aml+ball_female", "aml+ball_male", "tall+b_female","tall+b_male")
levels_target <- c("aml+nbl+ball_Female", "aml+nbl+ball_Male", "os+tall_Female","os+tall_Male")
load(paste(wd,"RNA/data/annotLookup.RData", sep="/"))

#INPUT: threshold -> float: The fdr cutoff value to use 
#OUTPUT: list containing the gene sets and their cat scores of 5 models
#       [[1]] -> Scores of the PMC SDA model trained on cancer types
#       [[c(2,3)]] -> PMC and TARGET group model results 
#       [[c(4,5)]] ->> PMC and TARGET sex-combined model results
read_output_files <- function(threshold = 0.8) {
  #Read in the files created by the pipeline containg the selected genes per model and their fdr value
  #first column are the gene names (NCBI ID)
  #Than the cat score for each group per gene is stated
  #lfdr < 0.8 for FNDR 
  #lfdr < 0.2 for FDR 
  #lfdr < thershold for any other desirable cutoff
  #last columns define if the gene is up or downregulated compared to the pooled mean. This is usefull since 
  #The decorelated nature of the CAT score, so a negative cat score doesn't mean that the gene is downregulated 
  #although this is true in most cases not for all 
  #The paper shows how we can obtain fold change back from cat-score, this was performed in the pipeline
  #here only the direciton (+1 or -1) for up or downregulation is given,
  #If you really want to see a specific effect for the gene, plot the boxplot from the normalised count data
  

  #Read output files 
  genes_group_gender_pmc <- read.csv("RNA/SDA/output/PMC_sign_genes_group_gender_all", col.names = 
                                       c("Gene", "index", "score", "cat_c_f", "cat_c_m", "cat_d_f", "cat_d_m", "lfdr", "HC", "d_c_f", "d_c_m", "d_d_f", "d_d_m"))
  genes_group_gender_target <- read.csv("RNA/SDA/output/TARGET_sign_genes_group_gender_all", col.names=
                                          c("Gene", "index", "score", "cat_c_f", "cat_c_m", "cat_d_f", "cat_d_m", "lfdr", "HC", "d_c_f", "d_c_m", "d_d_f", "d_d_m"))
  
  #remove genes below threshold
  genes_group_gender_pmc <- genes_group_gender_pmc[genes_group_gender_pmc$lfdr < threshold, ]
  genes_group_gender_target <- genes_group_gender_target[genes_group_gender_target$lfdr < threshold, ]
  
  genes_group_pmc <- read.csv("RNA/SDA/output/PMC_sign_genes_group_all", col.names =c("Gene", "index", "score", "cat_c", "cat_d", "lfdr", "HC", "d_c", "d_d"))
  genes_group_target <- read.csv("RNA/SDA/output/TARGET_sign_genes_group_all", col.names=c("Gene", "index", "score", "cat_c", "cat_d", "lfdr", "HC", "d_c", "d_d"))
  
  genes_disease_pmc <- read.csv("RNA/SDA/output/PMC_sign_genes_disease", col.names =c("Gene", "index", "score", "cat_AML", "cat_bALL","cat_burkitt", "cat_tALL", "lfdr", "HC", "d_AML","d_bALL","d_Burkitt","d_tALL"))
  
  #Remove genes below threshold
  genes_group_pmc <- genes_group_pmc[genes_group_pmc$lfdr < threshold, ]
  genes_group_target <- genes_group_target[genes_group_target$lfdr < threshold, ]
  
  genes_disease_pmc <- genes_disease_pmc[genes_disease_pmc$lfdr < threshold,]
  
  #Get rid of the extension to the gene id
  genes_group_gender_pmc$Gene <- substr(genes_group_gender_pmc$Gene,1,15)
  genes_group_pmc$Gene <- substr(genes_group_pmc$Gene,1,15)
  genes_disease_pmc$Gene <- substr(genes_disease_pmc$Gene,1,15)
  
  
  #Print the number of genes below the threshold of the group models
  nGenes_pmc <- dim(genes_group_pmc)[1]
  nGenes_target <- dim(genes_group_target)[1]
  
  print(paste("Number of genes pmc group:" , nGenes_pmc))
  print(paste("Number of genes target group:" , nGenes_target))

  #Print the number of genes below the threshold of the sex-combined models

  nGenes_pmc <- dim(genes_group_gender_pmc)[1]
  nGenes_target <- dim(genes_group_gender_target)[1]
  
  print(paste("Number of genes pmc group_gender:" , nGenes_pmc))
  print(paste("Number of genes target group_gender:" , nGenes_target))
  
  
  return(list("PMC_disease" = genes_disease_pmc,
              "PMC_group" = genes_group_pmc, 
              "TARGET_group" = genes_group_target, 
              "PMC_group_gender" = genes_group_gender_pmc,
              "TARGET_group_gender" = genes_group_gender_target
              ))
}  

#INPUT: threshold -> float: The fdr cutoff value to use 
#OUTPUT: list containing the gene sets and their cat scores of 5 models
#       [[c(1,2)]] -> PMC and TARGET cat scores
#       [[c(4,5)]] ->> PMC and TARGET sorted cat scores for difference between low and high irr boys 
important_genes_tables <- function(threshold=0.8) {
  data <- read_output_files(threshold)
  genes_group_pmc <- data[["PMC_group"]]
  genes_group_target <- data[["TARGET_group"]]
  
  genes_group_gender_pmc <- data[["PMC_group_gender"]]
  genes_group_gender_target <- data[["TARGET_group_gender"]]
  
  #Check the total overlap of the pmc and target group model
  pmc_target_genes_group <- genes_group_pmc %>% merge(genes_group_target, by="Gene")
  pmc_target_genes_group <- pmc_target_genes_group %>% merge(annotLookup, by.x="Gene", by.y="ensembl_gene_id", all.x=T)
  
  
  #Code that sorts the PMC model genes to where the low and high IRR boys differ the most from each other
  PMC_sex_combined <- (genes_group_gender_pmc %>%
         merge(annotLookup, by.x="Gene", by.y="ensembl_gene_id", all.x=T) %>%
          #Cat c_m is the cat score of the control (low IRR) male
          # cat_d_m is are the cat scores for high IRR males
         mutate(
           directed_cat_c_m = abs(cat_c_m) * d_c_m,
           directed_cat_d_m = abs(cat_d_m) * d_d_m,

           #Below the sort is just where low_IRR boys differ the most from high_IRR boys
           sort = (directed_cat_c_m - directed_cat_d_m)^2
         ))

  TARGET_sex_combined <- (genes_group_gender_target %>%
         merge(annotLookup, by.x="Gene", by.y="ensembl_gene_id", all.x=T) %>%
         mutate(
           directed_cat_c_m = abs(cat_c_m) * d_c_m,
           directed_cat_d_m = abs(cat_d_m) * d_d_m,


           #Below the sort is just where low_IRR boys differ the most from high_IRR boys
           sort = (directed_cat_c_m - directed_cat_d_m)^2

         )

         )


 return(list(
    "PMC_Group" = genes_group_pmc %>%  merge(annotLookup, by.x="Gene", by.y="ensembl_gene_id", all.x=T),
    "TARGET_Group" = genes_group_target %>% merge(annotLookup, by.x="Gene", by.y="ensembl_gene_id", all.x=T),
    "PMC_sex_combined" =  PMC_sex_combined,
    "TARGET_sex_combined" = TARGET_sex_combined
    ))
}


#Wrapper funciton for KEGG overrepresentation
#INPUT: genes -> vector: ENSEMBL ID's of genes with chosen FDR cutoff found by SDA
#       all_genes -> vector: ENSEMBL IDS of all genes in the study used as universe
#OUTPUT: Enrich KEGG result object and enrich plot is saved
over_representation <- function(genes, all_genes) {
  mapper <- bitr(genes,  "ENSEMBL", "ENTREZID", org.Hs.eg.db)
  universe <-  bitr(all_genes,  "ENSEMBL", "ENTREZID", org.Hs.eg.db)

  kegg_over<-enrichKEGG(
    mapper$ENTREZID, 
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe$ENTREZID,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.3,
    use_internal_data = FALSE
  )
  
  n <- sum(kegg_over@result$p.adjust < 0.05)
  
  print(enrichplot::dotplot(kegg_over, showCategory = n, font=9))
  
  return(kegg_over)
}

#Call this function to create over representation plots 
over_representation_analyses <- function(PMC_available=FALSE) {

  #Load the data, set threshold to 0.2
  data <- read_output_files(0.2)

  #Get the genes of the PMC group and sex-combined model
  if(PMC_available) {
    pmc_data <- load_PMC_data()
    all_genes_pmc <- substr(colnames(pmc_data[["X"]]),1,15)
    or_pmc_group<-over_representation((data[["PMC_group"]] %>% filter(d_c_m != d_d_m))$Gene, all_genes_pmc)
    or_pmc_sex_combined <- over_representation((data[["PMC_group_gender"]] %>% filter(d_c != d_d))$Gene, all_genes_pmc)
  }
 

  #TARGET group gender
  target_data <- load_TARGET_data()
  all_genes <- colnames(target_data[["X"]])
  over_representation((data[["TARGET_group"]]  %>% filter(d_c != d_d))$Gene, all_genes)
  tover_representation((data[["TARGET_group_gender"]] %>% filter(d_c_m != d_d_m))$Gene, all_genes)
}


#Download the kegg pathways
kegg <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")


#Implementation for PMC dataset, same can be performed for the TARGET dataset
kegg_interpretation_permutation <- function(dataset="PMC") {
  data <- read_output_files(0.2)
  
  if(dataset == "PMC") {
      #Do to parallelization of the HPC we splitted the 10000 permutations in 2 files. Load both and merge them
      permuted_kegg_scores.1 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/pmc_kegg_permutation_5001.csv", sep="/"))
      permuted_kegg_scores.2 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/pmc_kegg_permutation_5002.csv", sep="/"))

        #for all kegg pathways get the T2 value 
        t2_list <- tsquare_kegg(data[["PMC_group"]], "cat_d", "d_d")
  } else {
     permuted_kegg_scores.1 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/target_kegg_permutation_5001.csv", sep="/"))
     permuted_kegg_scores.2 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/target_kegg_permutation_5002.csv", sep="/"))

    #for all kegg pathways get the T2 value 
      t2_list <- tsquare_kegg(data[["TARGET_group"]], "cat_d", "d_d")
  }
  
  colnames(permuted_kegg_scores.1)[1] <- "ID"
  colnames(permuted_kegg_scores.2)[1] <- "ID"
  permuted_kegg_scores <- merge(permuted_kegg_scores.1, permuted_kegg_scores.2, by="ID")
  

  T2_values <- sapply(t2_list, function(a){as.double(a[["T2"]])})
  
  #Create dataframe with null pathway to have it filled already
  result <- data.frame("Id" = c("hsa000"), "name" = c("test"), "p_value"=c(1), "T2"=c(1), "p_adj"=c(1))
  
  #for all T2 values check significance
  for(name in names(T2_values)) {

    #Get the permuteted scores for this pathway, skip first column since its the name
    perm_distribution <- as.numeric(permuted_kegg_scores[permuted_kegg_scores[,1] == name,-1])
    
    #Calculate the percentage of scores from the permutations that is equal or greater than the found score by the real model for the pathway, than calculate the probability of finding our score 
    n_equal_or_higher_to_T2_value <- sum(perm_distribution >= T2_values[name])
    probability <- n_equal_or_higher_to_T2_value / length(perm_distribution)
    
    #Write the score of the pathway to the result dataframe, give p_adj 0 for now  
    result <- rbind(result, c("Id"= name,"name"= "", "p_value"=probability, "T2"=as.numeric(T2_values[name]), "p_adj"=0))
  }

  result <- result[-1,]
  #Add not only the ids but also the names of the KEGG pathway to the result dataframe
  result$name <- (kegg_pathway_id_to_name(result$Id))$to
  #Convert p values to doubles 
  result$p_value <- as.double(result$p_value)
  #Calculate adjusted p scores using Benjamini Hochberg method
  result$p_adj <- p.adjust(result$p_value, method="BH")

  result$T2 <- as.double(result$T2)
  
  #Return the dataframe
  return(result) 
}

#Function used to calcualte all the T^2 values given the found cat scores by SDA
tsquare_kegg <- function(df, group, dir="") {
  #some genes are then duplicated (only 1 or 2) so only keep the first one
  if(length(which(duplicated(df$Gene))) > 0) {
    df <- df[-which(duplicated(df$Gene)),]
  } 

  #Map ensembl ids to entrezids, as used by KEGG
  mapper <- bitr(df$Gene, "ENSEMBL", "ENTREZID", org.Hs.eg.db)

  list_T2<- list()
  #for every kegg pathway...
  for(pathway in unique(kegg$KEGGPATHID2EXTID$from)) { 
    #Get the entrez ids of the pathway
    entrez_ids_of_pathway <- kegg$KEGGPATHID2EXTID[kegg$KEGGPATHID2EXTI$from == pathway,"to"]
    ensembl_ids <- mapper[mapper$ENTREZID %in% entrez_ids_of_pathway, "ENSEMBL"]
    ensembl_ids <- ensembl_ids[!duplicated(ensembl_ids)]
    
    #Caluculate T^2 as descrbed by Zuber, V. and Strimmer, K. (2009).
    T2 <- sum(df[df$Gene %in% ensembl_ids, group]^2 )

    # As described we van also use the gene specific score for classification, calculate it, might be used in the future
    gene_specific_score <- data.frame("gene" = ensembl_ids, 
                                      "value" = sign(df[df$Gene %in% ensembl_ids, group]) * sqrt(T2)
                                      )
    #Calculate directed scores if you want to; Than give the name of the column that holds the value for transforming the cat scores to directed cat scores 
    directed_scores <- c()
    if(dir != "") {
      genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group, dir)]
      directed_scores <- abs(genes_of_intrest[,group]) * genes_of_intrest[, dir]
      names(directed_scores) <- genes_of_intrest$Gene
    }
    #Save all the genes that where used in the analysis and store there gene specific score with it 
    genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group)]
    
    list_T2[[pathway]] <- list("T2"= T2, "signed_score"=gene_specific_score, "cat_scores_group"= directed_scores)
  }
  
  return(list_T2)
}

kegg_pathway_id_to_name <- function(ids) {
  return(kegg$KEGGPATHID2NAME %>% dplyr::filter(from %in% ids))
}


#Helper function to order a list containg T^2 values
order_list <- function(ls){ 
  return(ls[order(sapply(ls, function(l) {unlist(l[1])}, simplify=T ), decreasing=T )])
}

kegg_interpretation_permutation("TARGET")


#Function to view a kegg pathway
kegg_pathway_viewer <- function(pathway) {
  data <- read_output_files(0.8)
  pmc_group <- data[["PMC_group"]]
  target_group <- data[["TARGET_group"]]
  
  pmc_group_gender <- data[["PMC_group_gender"]]
  target_group_gender <- data[["TARGET_group_gender"]]
  
  pmc_list_all <- tsquare_kegg(pmc_group, "cat_d", "d_d")
  ordered_pmc_all <-  order_list(pmc_list_all)
  
  genes<-pmc_list_all[[pathway]]$cat_scores_group
  idx <-genes
  
  mapper<-(bitr(names(idx), "ENSEMBL", c("ENTREZID", "UNIPROT" ), org.Hs.eg.db))
  names(idx)<- mapper[!duplicated(mapper$ENSEMBL),2]
  
  ss <- summary(idx)
  pathview(gene.data=idx,species="hsa",pathway.id=pathway, limit = c(-10, 10),  high = "#D40000FF", mid = "#FFFFFF", low = "#005EB8FF", cpd.idtype = "kegg",gene.idtype =
             "entrez", kegg.dir="kegg_analysis/kegg_data")
  
}


ls <- kegg_interpretation_permutation()
print(ls)


#Calculate P scores for the models from the 405 permuted 
calc_P <- function(df) {
  actual_kfold <- df$df_model_stats[1]
  actual_valid <- df$df_model_stats[2]
  print(actual_kfold)
  print(actual_valid)
  temp_kfold <- as.numeric(df[1,3:dim(df)[2]])
  temp_valid <- as.numeric(df[2,3:dim(df)[2]])

  p_kfold <- pnorm(actual_kfold, mean= mean(temp_kfold), sd=sd(temp_kfold), lower.tail = FALSE)
  p_valid <- pnorm(actual_valid, mean= mean(temp_valid), sd=sd(temp_valid), lower.tail=FALSE)
  
  return(list("kfold"= p_kfold, "valid"=p_valid))
}

#Example code using target group model, same holds for sex combined model and all  PMC models
target_gene_perm_label_group <- read.csv("RNA/Permutation_Tests/output/output_patient_permutation/group/TARGET_Group_model_stats_400.csv") 
target_model_perm_genes_group <- read.csv("RNA/Permutation_Tests/output/output_random_gene_permutation/TARGET_group_gene_perm_model_stats_400.csv")

calc_P(target_gene_perm_label_group)
calc_P(target_model_perm_genes_group)



#Example code on TARGET data, in paper PMC was used. Data can be obtained upon request
sda_pathway_plot <- function() {
  target_data <- load_TARGET_data()
  load("/Users/stevenwijnen/surfdrive/PMC/project/WGS/EXIT/kegg_t2.RData")
  scores <- read_output_files(0.8)
  scores[["TARGET_group_gender"]][,c(5,7)]
  t_genes <- scores[["TARGET_group_gender"]][scores[["TARGET_group_gender"]]$Gene %in% all[["hsa05165"]],]
t_genes.2 <- target_data[["X"]][,t_genes$Gene]
x <- apply(t_genes.2, 1, function(x) sum(x * t_genes[,c(5)]))
y <- apply(t_genes.2, 1, function(x) sum(x * t_genes[,c(6)]))

df <- data.frame(x,y)
df <- df %>% merge(target_data[["Y"]], by.x="row.names", by.y="patient_id")

ggplot(df , aes(x=x,y=y, color=Group_Gender, label=Row.names)) + geom_point() + 
  labs(x = "Component 1", y="Component 2") + coord_cartesian(xlim=c(-6000,10), ylim=c(-10,4000))+
 annotate("text", label="hsa05165", x=-250000,y=1100000, color="black") +
  scale_color_discrete(name="Class",labels=c('Low IRR female', 'Low IRR male', 'High IRR female', "High IRR male"))
}

sda_pathway_plot()