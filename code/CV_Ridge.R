# Function to implement Cross Validation for Ridge estimate
# Y nxs response matrix
# X nxp predictor matrix 
# lam vector of positive tuning parameters, default grid is 2^seq(-3, 3, 1).
# kfold is the number of folds for cross validation.
# The returned list includes
# lam--> optimal tuning parameter.
# min.error--> minimum average cross validation error for optimal parameters.
# avg.error--> average cross validation error across all folds.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple ridge
multi.ridge.cv <-function(X,Y, lambdaB, kfold=3, silent=TRUE){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]
  # Setting the tuning parameter values for RR
  lambdaB = sort(lambdaB)
  #This vector (CV_errors) will contain the prediction errors for each value of tuning parameter
  #lambdaOmega on the rows and lambdaB on the columns, when fitting the model.
  
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  #opt = list()
  opt.i=c()
  opt.lamB.out = c()
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]   # training predictors
    Y.tr.out = Y[i.tr.out,]    # training responses
    
    X.val.out = X[i.val.out,] # validation predictors
    Y.val.out = Y[i.val.out,]  # validation responses
    
    #Inner loop for hyper params  
    n.iner = nrow(Y.tr.out)
    K.iner = 3
    d.iner = ceiling(n.iner/K.iner)
    set.seed(123)
    i.mix.iner = sample(1:n.iner)
    folds.iner = vector(mode="list", length=K.iner) 
    for (k in 1:(K.iner-1)) {
      folds.iner[[k]] = i.mix.iner[((k-1)*d.iner+1):(k*d.iner)]
    }
    folds.iner[[K.iner]] = i.mix.iner[((K.iner-1)*d.iner+1):n.iner]  
    
    #This matrix (CV_errors) will contain the prediction errors for each combination of tuning parameters
    #lambdaOmega on the rows and lambdaB on the columns, when fitting the model.
    CV_errors.iner = rep(NA, length(lambdaB))
    CV_errors.val.iner = rep(list(CV_errors.iner), K.iner)

    for (k in 1:K.iner) {
      cat("inner Fold",k,"\n")
      
      i.tr.iner = unlist(folds.iner[-k])
      i.val.iner = folds.iner[[k]]
      
      X.tr.iner  = X.tr.out[i.tr.iner,]   # training predictors in the inner loop 
      Y.tr.iner  = Y.tr.out[i.tr.iner,]    # training responses in the inner loop 
      X.val.iner = X.tr.out[i.val.iner,] # validation predictors in the inner loop 
      Y.val.iner = Y.tr.out[i.val.iner,]  # validation responses in the inner loop 
      
      # Next, use the function JointEstim on the training inner data to get the 
      # L21-joint regression coefficients for all pairs of tuning parameter values in
      # lambdaO and lambdaB.
      # loop over all tuning parameters
      
        for(i in 1:length(lambdaB)){
          # compute the penalized regression matrix estimate with the training set
          B_hat.tr.iner = multi.ridge(X.tr.iner, Y.tr.iner, lambdaB[i]) 
          CV_errors.val.iner[[k]][i] = mean(mse.matrix(X.val.iner%*%B_hat.tr.iner,Y.val.iner))#eventually use the RV coef
        }
      
      idx.min=which.min((sapply(CV_errors.val.iner, min)))
      opt = which.min(CV_errors.val.iner[[idx.min]])      
      opt.i[k] = opt
      opt.lamB.out[k] = lambdaB[opt.i[k]]
      
    }
    
    #Now with the optimal tuning parameters from the inner loop
    #we validate on the left aside fold in the outer loop
    # compute the penalized regression matrix estimate with the training set
        B_hat.tr.out = multi.ridge(X.tr.out, Y.tr.out, opt.lamB.out[l]) 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%B_hat.tr.out, Y.val.out, na.rm = TRUE)#needed to compute the average metrics over all val folds
        CV_corr.val.out[[l]] = cor(X.val.out%*%B_hat.tr.out, Y.val.out)
  }
  
  best.idx.out=which.min((sapply(CV_errors.val.out, min)))
  best.lamB.final = opt.lamB.out[best.idx.out]
  #output on all data
  B_hat_final = multi.ridge(X, Y, best.lamB.final) 
  # return best lambdaB and other meaningful values
  return(list(best.lamB.final=best.lamB.final, B_hat_final=B_hat_final, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
