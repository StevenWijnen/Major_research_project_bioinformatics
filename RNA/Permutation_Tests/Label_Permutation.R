#In this script the labels are permuted 405 times;
#Foreach permutation a 5fold CV and validation score is calculated 
#The permutation scores are put together in the dataframe which is than saved as csv file and used 
#In further analysis

library(caret)
library(sda)
library(utils)
library(foreach)
library(doParallel)
library(pryr)
library(crossval)
library(doRNG)

wd <- getwd()
source(paste(wd, "RNA/Scripts/LoadData.R", sep="/"))
source(paste(wd, "RNA/Scripts/Crossval.R", sep="/"))


#INPUT: dataset -> string: either PMC or TARGET dataset to use form load_TARGET_data() or load_PMC_data()
#     : model -> string: Which model (Group or Group_Gender)
#     : dataset -> string: helper for the save funciton to add the name to the file, can be omitted
#     : levels -> vector: select the columns you want to use; Can be omitted 
#     : output_dir->  string: path to save to output
#OUTPUT: Writes the 400 permutations to a csv file with the fitrst entry the real model.        
main <- function(dataset, model, levels=c(), output_dir="") {
 set.seed(1816) 
  if(dataset == "PMC") {
    print("starting PMC")
    data <- load_PMC_data()
  }
  else {
    print("starting TARGET")
    data <- load_TARGET_data()
  }

 X <- data[["X"]]
 y <-  unlist(data[["Y"]][,model])
  
 y_model <- factor(y, levels=levels)
 inTrain <- createDataPartition(y_model, p = 0.75, list = FALSE)[,1]
  
  #train and test only for real model
  Xtrain <- as.matrix(X[ inTrain, ])
  Xvalid  <- as.matrix(X[-inTrain, ])
  
  ytrain <- (y_model[ inTrain])
  yvalid  <- (y_model[-inTrain])
  
  #Save real model scores
  df_model_stats <- model_pipeline(Xtrain, ytrain, Xvalid, yvalid, "real_model")
  
 print(paste("Number of cores ", detectCores()))
  
 finaldf <- parallel_func(X, y_model, dataset)
 print(finaldf) 
 print(df_model_stats) 
 write.csv(cbind(df_model_stats, finaldf), file=paste(output_dir, paste(dataset,"model_stats_400.csv", sep="_"), sep=""))

return(0)  
}

#Parallel computing is required to speed up the process
parallel_func <- function(X, y_model, dataset) {
  print(paste("Memory in use", mem_used()))
  start_time <- Sys.time()
  #Create a cluster, TODO adjust deteccores()-1 if you want to use more or less cores
   cl <- makeCluster(detectCores()-1, methods = FALSE, outfile="./log.txt", type="FORK")
  
  #The functions that need to be available on each thread 
  clusterExport(cl, c("model_pipeline", "crossVal", "predfun"))
  registerDoParallel(cl)
  
  #This is needed for randomization in parallel
  registerDoRNG(seed = 1616)
		
  #Create a temp df, a finite sized matrix speeds up computation
  temp_df <- matrix(data = NA, nrow = 2, ncol = 405)
  temp_df <- foreach(i=1:405, .combine=cbind, .errorhandling="remove", .packages = c("caret", "sda", "crossval", "pryr")) %dopar% {
    print("starting.....")
    #SBreak the correlation between the labels and the features
    yhat <- sample(y_model)
    #Create a new test valid split for each round
    parInTrain <- createDataPartition(yhat, p = 0.75, list = FALSE)[,1]
    
    parXtrain <- as.matrix(X[ parInTrain, ])
    parXvalid  <- as.matrix(X[-parInTrain, ])
    
    yhat_train <- (yhat[ parInTrain])
    yhat_valid  <- (yhat[-parInTrain])

    #Run the pipeline for the random model and retrieve kappa scores
    tempMatrix <-  model_pipeline(parXtrain, yhat_train, parXvalid, yhat_valid, i, dataset)
    

    tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  #Stop the cluster and report end time
  stopCluster(cl)
  end_time <- Sys.time()
  
  print(paste("time spent...:", (end_time - start_time)))

  return(temp_df)
}

#The pipeline to perform 5fold cv train the model on the train set and calculate validation set kappa statistic
model_pipeline <- function(Xtrain, ytrain, Xvalid, yvalid, run, dataset="") {
  #calculate ranking
  ra <- sda.ranking(Xtrain, ytrain, diagonal=FALSE, plot.fdr = FALSE, ranking.score = "entropy", verbose=FALSE)

  #Important features FNDR cutoff  
  nFeatures <- sum(ra[,"lfdr"] < 0.8)
  selVars = ra[,"idx"][1:nFeatures]
  
  #kFold CV to calculate accuracy
  model_stats <- crossVal(Xtrain, ytrain, nFeatures, K=5, B=10)
  
  sda.out<- sda(Xtrain[, selVars], ytrain, diagonal=FALSE, verbose=FALSE)
      
  pred_class <- predict(sda.out,  Xvalid[, selVars])$class
  cm <-caret::confusionMatrix(pred_class, yvalid)
  
  kappa_score <-  cm$overall[2]


  model_stats <- c("kFold"=model_stats, "test_set"=kappa_score)
 
  return(model_stats)
}


#Runs for different models:
#TODO umcomment to run model

#Group models
#main("PMC","Group", "PMC_Group", c("aml+ball", "tall+b"),"")
#main("TARGET", "Group", "TARGET_Group", c("aml+nbl+ball", "os+tall"), "")

#Sec-combined models
#main("PMC", "Group_Gender", "PMC_Group_Gender", c("aml+ball_female", "aml+ball_male", "tall+b_female","tall+b_male"),"")
main("TARGET", "Group_Gender", "TARGET_Group_Gender", c("aml+nbl+ball_Female", "aml+nbl+ball_Male", "os+tall_Female","os+tall_Male"))


