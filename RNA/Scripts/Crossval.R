library(crossval)
library(caret)

crossVal <- function(X,Y, nFeatures, K=5, B=10) {
  cv.sda = crossval(predfun, X, Y, K=K, B=B, numVars=nFeatures, verbose=FALSE)
  
  kappa_score <- cv.sda$stat
  
  return(kappa_score)
}

#Crossvalidation function needs working on
predfun <- function(Xtrain, Ytrain, Xtest, Ytest, numVars)
{
  # estimate ranking and determine the best numVars variables
  ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, verbose=FALSE, diagonal=FALSE)
  selVars = ra[,"idx"][1:numVars]
  
  sda.out = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=FALSE, verbose=FALSE)
  ynew = predict(sda.out, Xtest[, selVars, drop=FALSE], verbose=FALSE)$class
  
  cm <-caret::confusionMatrix(ynew, Ytest)
  kappa <- cm$overall[[2]]
  
  average_kappa <- mean(kappa)
  
  return(average_kappa)
}
