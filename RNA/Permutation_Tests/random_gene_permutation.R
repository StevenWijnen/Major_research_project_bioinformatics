#For this script we calculate the most important genes (using a training set, sda and FDR)
#Than we create an accuaracy distribution using kFold cross validation. Using this set of genes

library(sda)
library(caret)
library(crossval)
library(foreach)

wd <- getwd()
source(paste(wd, "RNA/Scripts/LoadData.R", sep="/"))


#INPUT: dataset -> string: either PMC or TARGET 
#     : model -> string: Group or Group_Gender 
#     : levels -> vector: LEvels of the classification
#     : output_dir -> string: path to output directory
#OUTPUT: saves the values for CV kappa and validation set Kappa to output file for the real model and all permuted models.
main <- function(dataset, model ,levels=c(), output_dir="") {
  set.seed(199816)
  
  #Load data
  if(dataset == "PMC") {
    print("starting PMC")
    data <- load_PMC_data()
  }
  else {
    print("starting TARGET")
    data <- load_TARGET_data()
  }


  X <- data[["X"]]
  Y <- data[["Y"]]

  #Get the right column (either Group or Group_Gender)
  y_model <- unlist(Y[,model])
  y_model <- factor(y_model, levels=levels)
  inTrain <- createDataPartition(y_model, p = 0.75, list = FALSE)[,1]
  
  #Create data parition
  Xtrain <- as.matrix(X[ inTrain, ])
  Xvalid  <- as.matrix(X[-inTrain, ])
  
  ytrain <- (y_model[ inTrain])
  yvalid  <- (y_model[-inTrain])
  
  #Train the real model
  df_model_stats <- model_pipeline(Xtrain, ytrain, Xvalid, yvalid, TRUE)

   #here we use, incontrary to the full model, the FDR to reduce the number of genes the random model can use
   #if we use the entire set, the probability of the random model to select a gene of intereset is high
   ra_2 <- sda.ranking(X, y_model, ranking.score="entropy", diagonal=FALSE)
   nGenes_d <- sum(ra_2[,"lfdr"] < 0.2)

  acc_random_model <- permute(Xtrain,ytrain, Xvalid, yvalid, nGenes_d)
  

  #Write to output file
  write.csv(cbind(df_model_stats, acc_random_model), 
            file= paste(paste(wd, output_dir, sep="/"), 
                             paste(dataset,"group_gene_perm_model_stats_400.csv", sep="_"), sep="/"))
  return(0)
}

#Input:  (Xtrain, ytrain, Xtest, ytest) -> dataframes: train test split data
#      : real_model-> boolean: do we need to choose random genes or do we use SDA to select genes 
#      : nGenes_d -> int: number of genes to pick for the random model
#Output 5CV and validation kappa scores
model_pipeline <- function(Xtrain, ytrain, Xvalid, yvalid, real_model=FALSE, nGenes_d=0) {
  ra <- sda.ranking(Xtrain, ytrain, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy", verbose=FALSE)
  if(real_model) {
    #For the real model use the correct genes
    nGenes <- sum(ra[,"lfdr"] < 0.8)
    sign_genes <- ra[,"idx"][1:nGenes]
  } else {
    #If the pipeline is called for a random model, sample random n genes
    #Supply the selected genes, either random or picked by SDA
    sign_genes <- ra[,"idx"][c(sample( 1:(dim(Xtrain))[2], nGenes_d))]
  }
  
  #Cross validation test score 
  cv.sda <- crossval(predfun, Xtrain, ytrain, genes=sign_genes, K=5, B=10, verbose=FALSE)

  #train the model for validation score 
  sda.out<- sda(Xtrain[, sign_genes], ytrain, diagonal=FALSE, verbose=FALSE)
  
  #Validaiton kappa score
  pred_class <- predict(sda.out,  Xvalid[, sign_genes])$class
  cm <-caret::confusionMatrix(pred_class, yvalid)
  
  kappa_score <-  cm$overall[2]

  model_stats <- c("kFold"=  cv.sda$stat, "test_set"=kappa_score)
  
  return(model_stats)
}

#predictor function used for the Crossvalidation implementation
#INPUT:   (Xtrain, ytrain, Xtest, ytest) -> dataframes: train test split data
#     : genes -> vector: containg the genes to use for SDA
predfun = function(Xtrain, ytrain, Xtest, ytest, genes = c())
{
  #No shrinkage, use all genes provided, since we want the model to use all genes we provided it with to make a corerct comparission
  #Here the shrinkage estimators are set really low such that all genes are used
  sda.out = sda(Xtrain[,genes], ytrain, diagonal=FALSE, verbose=FALSE, lambda=0.00001, lambda.var=0.00001, lambda.freqs=0.00001)
  
  #Used in the corss validation function, use the model trained on the ranom selected genes
  #and make a prediction on the test set
  ynew = predict(sda.out, Xtest[,genes], verbose=FALSE)$class
  
  cm <-caret::confusionMatrix(ynew, ytest)
  kappa <- cm$overall[[2]]
  
  #Since we perform 10 repetitions we pick the mean kappa score for the cross validation
  average_kappa <- mean(kappa)
  
  return(average_kappa)
}

#create a dataframe that keeps track of each itteration of a random model
#Here temp_df is filled with CV and validation kappa scores of the ranom models
permute <- function(Xtrain,ytrain, Xvalid, yvalid, genes) {
  temp_df <- matrix(data = NA, nrow = 2, ncol = 405)
  temp_df<-foreach(i= 1:405, .combine=cbind) %do% {
    print(i)
    temp_matrix <- model_pipeline(Xtrain,ytrain,Xvalid,yvalid, real_model=FALSE, nGenes_d = genes)

    temp_matrix
  }
  
  return(temp_df)
}

#Call the function

#main("PMC", "Group",  levels=c("aml+ball", "tall+b"), "RNA/Permutation_Tests/output")
#main("PMC","Group_Gender", levels= c("aml+ball_female", "aml+ball_male", "tall+b_female","tall+b_male"), "RNA/Permutation_Tests/output")
main("TARGET", "Group", levels = c("aml+nbl+ball", "os+tall"), "RNA/Permutation_Tests/output_test")
#main("TARGET","Group_Gender", levels = c("aml+nbl+ball_Female","aml+nbl+ball_Male", "os+tall_Female", "os+tall_Male"), "RNA/Permutation_Tests/output")
