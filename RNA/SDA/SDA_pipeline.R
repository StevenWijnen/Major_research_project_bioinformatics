#Script to run the SDA pipeline as described in the paper
#Here we ran the pipeline for TARGET and PMC datasets
library(sda)
library(caret)
library(crossval)

wd <- getwd()
#Load scripts
source(paste(wd, "RNA/Scripts/LoadData.R", sep="/"))
source(paste(wd, "RNA/Scripts/Crossval.R", sep="/"))

#INPUT: dataset to run the pipeline for (PMC or TARGET)
#     : the output directory in which you want to save the kappa scores and most important genes
main <- function(dataset, output_dir) {
  set.seed(1998)
  if(dataset == "PMC") {
    print("starting PMC")
    data <- load_PMC_data()
  }
  else {
    print("starting TARGET")
    data <- load_TARGET_data()
  }

  Y <- data[["Y"]] #patient data
  Y$Disease <- as.factor(Y$Disease)
  Y$Group <- as.factor(Y$Group)
  Y$Group_Gender <- as.factor(Y$Group_Gender)
  
  X <- data[["X"]] # count data 
  #Create a train validation split
  inTrain <- createDataPartition(Y$Disease, p=0.75, list = FALSE)[,1]

#train and validation sets
  Xtrain <- X[inTrain,]
  Xvalid <- X[-inTrain,]
  ytrain <- Y[inTrain,]
  yvalid <- Y[-inTrain,]
  
  #Perform SDA on cancer type
  sda_cancer_types <- model_pipeline(Xtrain, ytrain$Disease, Xvalid, yvalid$Disease)
  
  #extract significant genes for cancer type
  sign_genes_cancer_type <- sda_cancer_types[[1]][,1]
  #Save significant genes and kappa scores
  write.csv(sda_cancer_types[[2]], file= paste(output_dir, paste(dataset,"cancer_types_model_stats.csv", sep="_"), sep=""))
  write.csv(sda_cancer_types[[1]], file= paste(output_dir, paste(dataset,"sign_genes_disease", sep="_"),sep=""))
  
  #Code will continue with and without the cancer_type_genes (See scheme)
  #all genes; Group model
  sda_disease_vs_healthy_all_genes <- model_pipeline(Xtrain, ytrain$Group, Xvalid, yvalid$Group)
  
  #skip important genes form cancer type prediction; if this still gives a good accuracy than we know that 
  #the SDA disease_healthy model isn't just predicting cancer types
  sda_disease_vs_healthy_skip_genes <- model_pipeline(
    Xtrain[, -sign_genes_cancer_type], 
    ytrain$Group, Xvalid[,-sign_genes_cancer_type], 
    yvalid$Group
    )
  
  #Save kappa scores and cat scores
  write.csv(sda_disease_vs_healthy_all_genes[[2]], file=paste(output_dir, paste(dataset,"group_all_genes_model_stats.csv", sep="_"),sep="/"))
  write.csv(sda_disease_vs_healthy_skip_genes[[2]], file=paste(output_dir, paste(dataset,"group_subset_genes_model_stats.csv", sep="_"),sep="/"))
  
  write.csv(sda_disease_vs_healthy_all_genes[[1]], file= paste(output_dir, paste(dataset,"sign_genes_group_all", sep="_"), sep="/"))
  
  #Now on the sex combined model
  sda_group_gender <- model_pipeline(Xtrain, ytrain$Group_Gender, Xvalid, yvalid$Group_Gender)
  sda_group_gender_skip_genes <- model_pipeline(Xtrain[,-sign_genes_cancer_type], ytrain$Group_Gender, Xvalid[,-sign_genes_cancer_type], yvalid$Group_Gender)
  
  #save cat scores and kappa scores
  write.csv(sda_group_gender[[2]], file=paste(output_dir, paste(dataset,"group_gender_all_model_stats.csv", sep="_"),sep="/"))
  write.csv(sda_group_gender_skip_genes[[2]], file=paste(output_dir, paste(dataset,"group_gender_subset_genes_model_stats.csv", sep="_"),sep="/"))
  
  write.csv(sda_group_gender[[1]], file= paste(output_dir, paste(dataset,"sign_genes_group_gender_all", sep="_"),sep="/"))

  return(sda_group_gender[[1]]) 
}

#Pipeline run for each model
model_pipeline <- function(Xtrain, ytrain, Xvalid, yvalid) {
  print("starting pipeline....")
  #calculate ranking
  ra <- sda.ranking(Xtrain, ytrain, fdr=TRUE, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy")

  nFeatures <- sum(ra[,"lfdr"] < 0.8)
  
  #extract important features 
  selVars = ra[,"idx"][1:nFeatures]

  #kFold CV to calculate accuracy
  model_stats <- crossVal(Xtrain, ytrain, nFeatures, K=5, B=10)
  #train model for validation set
  sda.out <- sda(Xtrain[, selVars], ytrain, diagonal=FALSE)
  pred_class <- predict(sda.out,  Xvalid[, selVars])$class
  #report validation test score
  cm <-caret::confusionMatrix(pred_class, yvalid)
  kappa_score <-  cm$overall[2]

  #Save 5fold cv and validation kappa score
  model_stats <- c("kFold_kappa_score" = model_stats, "test_score" =  kappa_score)
  
  #For the final model we use the entire dataset to get the most accurate cat scores, no kappa scores are calculated
  #This step is for further analyses on the cat scores results
  tmp = centroids(rbind(Xtrain, Xvalid), c(ytrain, yvalid), var.groups = FALSE, centered.data = TRUE, verbose = FALSE)
  cl.count <- length(levels(ytrain))
  #Calculate the fold change, using pooled mean and gorup mean
  mu = tmp$means[, 1:cl.count]
  mup = tmp$means[, cl.count + 1]
  diff = sign(mu[,] - mup)
  
  #Get most important features
  ra <- sda.ranking(rbind(Xtrain, Xvalid), c(ytrain, yvalid), fdr=TRUE, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy")
  nFeatures <- sum(ra[,"lfdr"] < 0.8)

  ww <- ra[,][1:nFeatures,] %>% merge(diff, by = "row.names") %>% tibble::column_to_rownames(var="Row.names")
  
  #Return cat scores trained for the entire dataset and kappa scores obtained from train and validation set
  return(list(ww, model_stats))

}

#PMC is commented out since no data is publicly available 
#pmc_genes <- main("PMC")
target_genes <- main("TARGET", "/Users/stevenwijnen/surfdrive/PMC/Major_research_project_bioinformatics/RNA/SDA/output_test_small")
