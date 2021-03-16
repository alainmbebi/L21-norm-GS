#function for cross validation in Multitrait Bayes B using the R package BGLR
multi.BayesB.cv <-function(X,Y, kfold=5){
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

  #This vector (CV_errors) will contain the prediction errors for each value of tuning parameter
  #lambdaOmega on the rows and lambdaB on the columns, when fitting the model.
  SVD=list()
  U=list()
  D=list()
  V=list()
  BETA_MBayesB=list()
  Bhat_MBayesB=list()
  
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]   # training predictors
    Y.tr.out = Y[i.tr.out,]    # training responses
    
    X.val.out = X[i.val.out,] # validation predictors
    Y.val.out = Y[i.val.out,]  # validation responses
    
    # compute the penalized regression matrix estimate with the training set
    
    # Singular value decompositions Y=UDV' Bhat_MBayesB
    SVD[[l]]=svd(Y.tr.out)
    U[[l]]=SVD[[l]]$u
    D[[l]]=diag(SVD[[l]]$d)
    V[[l]]=SVD[[l]]$v
  
  
    BETA_MBayesB[[l]]=matrix(nrow=ncol( X.tr.out),ncol=ncol(Y.tr.out))
  
   ETA=list(list(X=X.tr.out, model='BayesB'))
  for(i in 1:ncol(Y.tr.out)){
    fm=BGLR(y=U[[l]][,i],ETA=ETA,verbose=F) #use more iterations!
    BETA_MBayesB[[l]][,i]=fm$ETA[[1]]$b
  }
  
  # Rotating coefficients to put them in marker space
  Bhat_MBayesB[[l]]=BETA_MBayesB[[l]]%*%D[[l]]%*%t(SVD[[l]]$v)
  
 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%Bhat_MBayesB[[l]], Y.val.out, na.rm = TRUE)#needed to compute the average metrics over all val folds
        CV_corr.val.out[[l]] = cor(X.val.out%*%Bhat_MBayesB[[l]], Y.val.out)
  }

  #output on all data
  # return best lambdaB and other meaningful values
  return(list(B_hat_final=Bhat_MBayesB, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}
