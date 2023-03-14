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
  #If you really want to see a specific effect for the gene, plot the boxplot from the normalised count data see: plot_boxplot() function
  

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


over_representation_analyses <- function(PMC_available=FALSE) {

  #Load the data, set threshold to 0.2
  data <- read_output_files(0.2)

  #Get the genes of the PMC group and sex-combined model
  if(PMC_available) {
    pmc_data <- load_PMC_data()
    all_genes_pmc <- substr(colnames(pmc_data[["X"]]),1,15)
    or_pmc_group<-over_representation(data[["PMC_group"]]$Gene, all_genes_pmc)
    or_pmc_sex_combined <- over_representation(data[["PMC_sex_combined"]]$Gene, all_genes_pmc)
  }
 

  #TARGET group gender
  target_data <- load_TARGET_data()
  all_genes <- colnames(target_data[["X"]])
  over_representation(data[["TARGET_group"]]$Gene, all_genes)

  over_representation(data[["TARGET_group_gender"]]$Gene, all_genes)
  
}


over_representation_analyses()





# kegg_pathway_id_to_name <- function(ids) {
#   return(kegg$KEGGPATHID2NAME %>% dplyr::filter(from %in% ids))
# }

# tsquare_kegg <- function(df, group, dir="") {
#   kegg <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")

#   #some genes are then duplicated (only 1 or 2) so only keep the first one
#   #apperently -c(0) returns empty dataframe, so first check if there are duplcites
#   if(length(which(duplicated(df$Gene))) > 0) {
#     df <- df[-which(duplicated(df$Gene)),]
#   } 
#   #Map ensembl ids to entrezids
#   mapper <- bitr(df$Gene, "ENSEMBL", "ENTREZID", org.Hs.eg.db)

#   list_T2<- list()
#   #for every kegg pathway...
#   for(pathway in unique(kegg$KEGGPATHID2EXTID$from)) { 
#     #Get the entrez ids of the pathway
#     entrez_ids_of_pathway <- kegg$KEGGPATHID2EXTID[kegg$KEGGPATHID2EXTI$from == pathway,"to"]
#     ensembl_ids <- mapper[mapper$ENTREZID %in% entrez_ids_of_pathway, "ENSEMBL"]
#     ensembl_ids <- ensembl_ids[!duplicated(ensembl_ids)]
#     T2 <- sum(df[df$Gene %in% ensembl_ids, group]^2 )

#         #make class which can hold specific genes if we need that 
#     gene_specific_score <- data.frame("gene" = ensembl_ids, 
#                                       "value" = sign(df[df$Gene %in% ensembl_ids, group]) * sqrt(T2)
#                                       )
#     directed_scores <- c()
#     if(dir != "") {
#       genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group, dir)]
#       directed_scores <- abs(genes_of_intrest[,group]) * genes_of_intrest[, dir]
#       names(directed_scores) <- genes_of_intrest$Gene
#     }
    
#     genes_of_intrest <- df[df$Gene %in% ensembl_ids, c("Gene", group)]
    
#     list_T2[[pathway]] <- list("T2"= T2, "signed_score"=gene_specific_score, "cat_scores_group"= directed_scores)
#   }
  
#   return(list_T2)
# }


# kegg_interpretation_sorted <- function(group1="PMC_group", group2="TARGET_group") {
#   data <- read_output_files(0.8)

#   list1 <-  tsquare_kegg(data[[group1]], "cat_d", "d_d")
#   ordered_list1 <-  order_list(list1)
#   list1_T2 <- sapply(ordered_list1, function(ls){return(ls$T2)})
  
#   list2 <- tsquare_kegg(data[[group2]], 'cat_d', "d_d")
#   ordered_list2 <- order_list(list2)
#   list2_T2 <- sapply(ordered_list2, function(ls){return(ls$T2)})
  
#   kegg_names1<- kegg_pathway_id_to_name(names(ordered_list1[1:20]))
#   df1 <- cbind(kegg_names1, "T2"=list1_T2[1:20])
  
#   kegg_names2 <-  kegg_pathway_id_to_name(names(ordered_list2[1:20]))
#   df2 <- cbind(kegg_names2, "T2"=list2_T2[1:20])
  
#   #TODO,Uncomment If you want to insepct dataframes in R studio
#   #View(df1)
#   #View(df2)
#   #View(merge(kegg_names1, kegg_names2, by="from"))
  
#   return(list("group1" = df1, "group2" = df2))
# }


# kegg_interpretation_permutation <- function() {
#   data <- read_output_files(0.8)
#   pmc_permuted_kegg_scores.1 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/kegg_scores_permutation/pmc_kegg_permutation_5001.csv", sep="/"))
#   pmc_permuted_kegg_scores.2 <- read_csv(paste(wd, "RNA/Permutation_Tests/output/output_kegg_permutation/kegg_scores_permutation/pmc_kegg_permutation_5002.csv", sep="/"))

#   colnames(pmc_permuted_kegg_scores.1)[1] <- "ID"
#   colnames(pmc_permuted_kegg_scores.2)[1] <- "ID"
#   pmc_permuted_kegg_scores <- merge(pmc_permuted_kegg_scores.1, pmc_permuted_kegg_scores.2, by="ID")
  
#   #first calculate T^2 values for pmc 
#   pmc_list <- tsquare_kegg(data[["PMC_group"]], "cat_d", "d_d")
#   #for all kegg pathways get the T2 value 
#   T2_values <- sapply(pmc_list, function(a){as.double(a[["T2"]])})
#   print(typeof(T2_values))
  
#   result <- data.frame("Id" = c("hsa000"), "name" = c("test"), "p_value"=c(1), "T2"=c(1), "p_adj"=c(1))
#   #for all T2 values check significance
#   for(name in names(T2_values)) {

#     perm_distribution <- as.numeric(pmc_permuted_kegg_scores[pmc_permuted_kegg_scores[,1] == name,-1])
#     n_equal_or_higher_to_T2_value <- sum(perm_distribution >= T2_values[name])
    
#     probability <- n_equal_or_higher_to_T2_value / length(perm_distribution)
    
#     result <- rbind(result, c("Id"= name,"name"= "", "p_value"=probability, "T2"=as.numeric(T2_values[name]), "p_adj"=0))
#   }

  
#   result <- result[-1,]
#   result$name <- (kegg_pathway_id_to_name(result$Id))$to
#   result$p_value <- as.double(result$p_value)
  
#   result$p_adj <- p.adjust(result$p_value, method="BH")

#   result$T2 <- as.double(result$T2)
  
#   return(result)
  
# }


# order_list <- function(ls){ 
#   return(ls[order(sapply(ls, function(l) {unlist(l[1])}, simplify=T ), decreasing=T )])
# }

# kegg_pathway_viewer <- function(pathway) {
#   data <- read_output_files(0.8)
#   pmc_group <- data[["PMC_group"]]
#   target_group <- data[["TARGET_group"]]
  
#   pmc_group_gender <- data[["PMC_group_gender"]]
#   target_group_gender <- data[["TARGET_group_gender"]]
  
#   pmc_list_all <- tsquare_kegg(pmc_group, "cat_d", "d_d")
#   ordered_pmc_all <-  order_list(pmc_list_all)
  
#   genes<-pmc_list_all[[pathway]]$cat_scores_group
#   idx <-genes
  
#   mapper<-(bitr(names(idx), "ENSEMBL", c("ENTREZID", "UNIPROT" ), org.Hs.eg.db))
#   names(idx)<- mapper[!duplicated(mapper$ENSEMBL),2]
  
#   ss <- summary(idx)
#   pathview(gene.data=idx,species="hsa",pathway.id=pathway, limit = c(-10, 10),  high = "#D40000FF", mid = "#FFFFFF", low = "#005EB8FF", cpd.idtype = "kegg",gene.idtype =
#              "entrez", kegg.dir="kegg_analysis/kegg_data")
  
# }



# data <- read_output_files(0.8)
# #PMC group and group gender
# all_genes <- substr(colnames(pmc_data[["X"]]),1,15)
# t1<-over_representation(data[["PMC_group"]]$Gene, all_genes)

# #TARGET group and group gender
# over_representation(data[["TARGET_group_gender"]]$Gene, colnames(target_data[["X"]]))




# #gseKEGG use the gene set enrichment analysis and enrichKEGG use over-representation test, 
# #and GESA doesn't need to run differential gene expression analysis beforehand

# genes_group_gender_target <- read.csv("Scripts/sda_pipeline_output/TARGET_sign_genes_group_gender_all", col.names=c("Gene", "index", "score", "cat_c_f", "cat_c_m", "cat_d_f", "cat_d_m"))

# sum(genes_group_gender_target$Gene %in% colnames(target_data[["X"]][,which(!colnames(target_data[["X"]])%in% substr(colnames(pmc_data[["X"]]), 1,15))]))



# tt <- pmc_data[["Y"]]
# X <- pmc_data[["X"]][tt$Biomaterial_Id,] 
# all(rownames(X) == tt$Biomaterial_Id)
# y <- tt$Group_Gender

# ra <- sda.ranking(X,y , fdr=TRUE, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy")
# ra.1 <- ra[ra[,"lfdr"] < 0.8, ]
# nFeatures <- sum(ra[,"lfdr"] < 0.8)

# selVars = ra[,"idx"][1:nFeatures]
# df <- rownames_to_column(as.data.frame(ra.1), var="Gene")
# df$Gene <- substr(df$Gene, 1,15)
# pmc_list <- tsquare_kegg(df, "cat.tall+b_male")
# ordered_pmc_males <- order_list(pmc_list)
# pmc_kegg_names_all <- kegg_pathway_id_to_name(names(ordered_pmc_males[1:20]))




# test <- tsquare_kegg(data[["PMC_group"]], "cat_d")
# group <- kegg_interpretation_sorted()[[1]]

# all <- list()
# for(name in names(test)) {
#   print(name)
#   if(name %in% group[1:20,"from"] | name %in% t1$ID) {
#     all[[name]] <- test[[name]]$signed_score$gene   
#   }
# }

# save(all, file="/Users/stevenwijnen/surfdrive/PMC/project/WGS/EXIT/kegg_t2.RData")

# library(rgl)
# test <- read_output_files(0.8)
# test[["PMC_group_gender"]][,c(5,7)]
# t_genes<-pmc_data[["X"]][,test[["PMC_group_gender"]]$Gene]
# x <- apply(t_genes, 1, function(x) sum(x * test[["PMC_group_gender"]][,c(4)]))
# y <- apply(t_genes, 1, function(x) sum(x * test[["PMC_group_gender"]][,c(6)]))
# z <- apply(t_genes, 1, function(x) sum(x * test[["PMC_group_gender"]][,c(6)]))

# df <- data.frame(x,y,z)
# df <- df %>% merge(pmc_data[["Y"]], by.x="row.names", by.y="Biomaterial_Id")

# ggplot(df , aes(x=x,y=y, color=Group_Gender, label=Row.names)) + geom_point() 


# plot_ly(x=df$x, y=df$y, z=df$z, type="scatter3d", mode="markers", color=df$Disease)

 
 
#  t_genes<-pmc_data[["X"]][,test[["PMC_group"]]$Gene]
#  x <- apply(t_genes, 1, function(x) sum(x * test[["PMC_group"]][,c(4)]))
#  y <- apply(t_genes, 1, function(x) sum(x * test[["PMC_group"]][,c(5)]))

#  df <- data.frame(x,y,z)
#  df <- df %>% merge(pmc_data[["Y"]], by.x="row.names", by.y="Biomaterial_Id")
#  ggplot(df, aes(x=x,y=y, color=Disease, label=Row.names)) + geom_point() 
 
#  X <- pmc_data[["X"]]
#  X <- X[which(!row.names(X) %in% weird_ones) ,]
#  Y <- pmc_data[["Y"]]
#  Y <- Y[which(!Y$Biomaterial_Id %in% weird_ones),]
#  ra <- sda.ranking(X,Y$Group_Gender)
#  ra.fndr <- ra[which(ra[,"lfdr"] < 0.8),]
#  t_genes<-pmc_data[["X"]][,row.names(ra.fndr)] %>% filter()
 
#  x <- apply(t_genes, 1, function(x) sum(x * ra.fndr[,4]))
#  y <- apply(t_genes, 1, function(x) sum(x * ra.fndr[,6]))

#  df <- data.frame(x,y,z)
#  df <- df %>% merge(pmc_data[["Y"]], by.x="row.names", by.y="Biomaterial_Id")
#  ggplot(df, aes(x=x,y=y, color=Disease, label=Row.names)) + geom_point()  

 
 
 
 
 