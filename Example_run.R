# '@Alain Mbebi
set.seed(123)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# '@ libraries
library(data.table)
library(MASS)
library(corpcor)
library(Matrix)           
library(lattice)
library(mvtnorm)
library(MRCE)             # This will help  to compare with the MRCE model of Rothman et.al 2010
library(remMap)           # This will help  to compare with the remMap model of Jie Peng et.al 2010
library(glasso)           # This will help to estimate precision matrix using Graphical Lasso
library(glassoFast)
library(matrixcalc)
library(CVTuningCov)
library(glmnet)
library(BGLR)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function for features selection algorithm 1
FeatSelect<-function(lambdaB, X, Y, eps=1e-6, zeta=1e-8){
  #--------------------------
  A=cbind(X, lambdaB*diag(nrow(Y)))
  #create an empty matrice to store the beta's coeficients to be estimated.
  B_l21fs=matrix(0, nrow=ncol(X), ncol=ncol(Y))
  
  #--------------------------
  # t=0 initialization
  diag_D_0=c(rep(1, (ncol(X)+nrow(Y))))                # Initialize D as the identity matrix of dim the number nrow of W_true
  invdiag_D_0=1/diag_D_0              		       # Digonal matrix representing the inverse of D_0
  
  #--------------------------
  #t=1
  W_1=round((invdiag_D_0*t(A))%*%solve(A%*%(invdiag_D_0*t(A)))%*%Y ,9)    
  l21W_1=c(rep(0, (ncol(X)+nrow(Y))))
  for(i in 1:(ncol(X)+nrow(Y))){
    l21W_1[i]=2*norm((W_1[i,]),type = "2")# + zeta
  }
  invdiag_D_1=l21W_1
  
  #--------------------------
  #Then the update of W at t=2 as:
  W_2=round((invdiag_D_1*t(A))%*%solve(A%*%(invdiag_D_1*t(A)))%*%Y ,9)   #that is t=2
  l21W_2=c(rep(0, (ncol(X)+nrow(Y))))
  for(i in 1:(ncol(X)+nrow(Y))){
    l21W_2[i]=2*norm((W_2[i,]),type = "2") #+ zeta
  }
  invdiag_D_2=l21W_2
  
  #--------------------------
  tcont=0
  while((sum(l21W_2)<=sum(l21W_1))==TRUE && (min(diag((A%*%(invdiag_D_0*t(A))))) > eps) && tcont<=itermax){
    W_1=W_2
    l21W_1=l21W_2                                                                    
    invdiag_D_1=invdiag_D_2                                                          
    
    W_2=round((invdiag_D_1*t(A))%*%solve(A%*%(invdiag_D_1*t(A)))%*%Y ,9)       
    l21W_2=c(rep(0,(ncol(X)+nrow(Y))))
    for(i in 1:(ncol(X)+nrow(Y))){
      l21W_2[i]=2*norm((W_2[i,]),type = "2") 
    }
    invdiag_D_2=l21W_2
    tcont <- sum(tcont, 1)
    
    print(tcont)
    B_l21fs=round(W_1[1:ncol(X),],9)                         		   
  }
  return(B_l21fs)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Cross Validation in feature selection algorithm 1
CV_feat_sel = function(Y, X, lam, kfold=5) {
  
  lam = sort(lam)
  # initialize
  CV_errors = array(0, c(length(lam), kfold))
  
  # designate folds and shuffle -- ensures randomized folds
  n = nrow(Y)
  ind = sample(n)
  
  # parse data into folds and perform CV
  
  for (k in 1:kfold) {
    
    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    
    # validation set
    Y.valid = Y[leave.out,,drop=FALSE ]
    X.valid = X[leave.out,,drop=FALSE ]
    
    #this centers the training and validation data, uncomment if needed
    #Y.train_bar = apply(Y.train, 2, mean)
    #Y.train = scale(Y.train, center = Y.train_bar, scale = FALSE)
    #Y.valid = scale(Y.valid, center = Y.train_bar, scale = FALSE)
    #Y.train = scale(Y.train, center = TRUE, scale = TRUE)
    #Y.valid = scale(Y.valid, attr(Y.train, "scaled:center"), attr(Y.train, "scaled:scale"))
    
    #X.train_bar = apply(X.train, 2, mean)
    #X.train = scale(X.train, center = X.train_bar, scale = FALSE)
    #X.valid = scale(X.valid, center = X.train_bar, scale = FALSE)
    
    #X.train = scale(X.train, center = TRUE, scale = TRUE)
    #X.valid = scale(X.valid, attr(X.train, "scaled:center"), attr(X.train, "scaled:scale"))

    
    # loop over all tuning parameters
    for (i in 1:length(lam)) {
      
      # compute the penalized regression matrix estimate with the training set
      B_hat.train=FeatSelect(lam[i], X.train, Y.train) 
      
      # compute the CV errors
      CV_errors[i,k] = CV_errors[i,k] + mean((Y.valid-X.valid%*%B_hat.train)^2)
    }
    
    
  }
  
  # determine optimal tuning parameters
  AVG = apply(CV_errors, 2, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  
  
  # return best lam and other meaningful values
  return(list(lam = best_lam, min.error = error, avg.error = AVG, cv.error = CV_errors))
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for multiple output regression (MOR)
cmor<-function(lambda1, lambda2, lambda3, lambda4, X, Y, eps=1e-6, zeta=1e-8, tol.out=1e-6){

# t=0 initialization
  
  Sigma_0=diag(ncol(Y))                           # Initialize Omega as the identity matrix of dim the number ncol of Y 
  Omega_0=diag(ncol(Y))                           # Initialize Sigma as the identity matrix of dim the number ncol of Y 
  invSigma0=solve(Sigma_0)
  invOmega0=solve(Omega_0)
  P0=t(chol(lambda1*Omega_0 +  lambda2*Sigma_0))  # to get the lower triangular matrix P (notice the transpose otherwise upper triangular)
  
  #----------------------
  SVD_P0=svd(P0)
  U2_0=SVD_P0$u
  Sigma2_0=diag(SVD_P0$d)
  Sigma2vect_0=c(SVD_P0$d)
  tV2_0=SVD_P0$v
  
  #----------------------
  SVD_X=svd(X)
  U1=SVD_X$u
  Sigma1=diag(SVD_X$d)
  Sigma1vect=c(SVD_X$d)
  tV1=t(SVD_X$v)
  
  #----------------------
  Btilde0= matrix(0, nrow=nrow(X), ncol=ncol(Y))
  S0=tV1%*%t(X)%*%Y%*%U2_0
  diagSigma1=(Sigma1vect)^2 #length p
  diagSigma2_0=(Sigma2vect_0)^2 #length s
  for(i in 1:(nrow(X))){
    for(j in 1:(ncol(Y))){
      Btilde0[i,j]<-S0[i,j]/(diagSigma1[i] + diagSigma2_0[j])
    }
  }
  B0=t(tV1)%*%Btilde0%*%t(U2_0)
  
  residual.cov_0=crossprod(Y-X%*%B0)/nrow(Y)
  obj_0=matrix.trace(residual.cov_0%*%invOmega0)-nrow(Y)*determinant(invOmega0, logarithm=TRUE)$mod[1] +
   lambda1*matrix.trace(t(B0)%*%B0) + lambda2*matrix.trace(invSigma0%*%t(B0)%*%B0) -ncol(X)*determinant(invSigma0, logarithm=TRUE)$mod[1] +
   lambda3*matrix.trace(invOmega0) +  lambda4*matrix.trace(invSigma0)
  #----------------------
  Sigma_1=(lambda2*t(B0)%*%B0 + lambda4*diag(ncol(Y)))/nrow(X)
  Omega_1=(t(Y-X%*%B0)%*%(Y-X%*%B0) + lambda3*diag(ncol(Y)))/nrow(Y)
  invSigma1=solve(Sigma_1)
  invOmega1=solve(Omega_1)
  P1=t(chol(lambda1*Omega_1 +  lambda2*Sigma_1)) 
  
  #----------------------
  SVD_P1=svd(P1)
  U2_1=SVD_P1$u
  Sigma2_1=diag(SVD_P1$d)
  Sigma2vect_1=c(SVD_P1$d)
  tV2_1=SVD_P1$v
  
  #----------------------
  Btilde1= matrix(0, nrow=nrow(X), ncol=ncol(Y))
  S1=tV1%*%t(X)%*%Y%*%U2_1
  diagSigma1=(Sigma1vect)^2 #length p
  diagSigma2_1=(Sigma2vect_1)^2 #length s
  for(i in 1:(nrow(X))){
    for(j in 1:(ncol(Y))){
      Btilde1[i,j]<-S1[i,j]/(diagSigma1[i] + diagSigma2_1[j])
    }
  }
  B1=t(tV1)%*%Btilde1%*%t(U2_1)
  
  residual.cov_1=crossprod(Y-X%*%B1)/nrow(Y)
  obj_1=matrix.trace(residual.cov_1%*%invOmega1)-nrow(Y)*determinant(invOmega1, logarithm=TRUE)$mod[1] +
   lambda1*matrix.trace(t(B1)%*%B1) + lambda2*matrix.trace(invSigma1%*%t(B1)%*%B1) -ncol(X)*determinant(invSigma1, logarithm=TRUE)$mod[1] +
   lambda3*matrix.trace(invOmega1) +  lambda4*matrix.trace(invSigma1)
  #----------------------
  
  tcontjmor=0
  conv_crit=tol.out*sum(diag(crossprod(Y))/nrow(Y))
  #while((sum(abs(B1))<=sum(abs(B0)))==TRUE  && tcontjmor<=itermax){ 
  while(((obj_0-obj_1)>conv_crit)==TRUE  && (min(diag(residual.cov_0)) > eps) && tcontjmor<=itermax){ 
    
    Sigma_0=Sigma_1
    Omega_0=Omega_1  
    invSigma0=invSigma1                                         
    invOmega0=invOmega1   
    P0=P1
    #----------------------
    #----------------------
    SVD_P0=SVD_P1
    U2_0=U2_1
    Sigma2vect_0=Sigma2vect_1
    tV2_0=tV2_1
    
    #----------------------
    Btilde0= Btilde1
    S0=S1
    B0=B1
        
    #----------------------
    residual.cov_0=residual.cov_1
    obj_0=obj_1
    
    #----------------------
    #----------------------
    
    SVD_P0=svd(P0)
    U2_0=SVD_P0$u
    Sigma2_0=diag(SVD_P0$d)
    Sigma2vect_0=c(SVD_P0$d)
    tV2_0=SVD_P0$v
    
    #---------------------- 
    Sigma_1=(lambda2*t(B0)%*%B0 + lambda4*diag(ncol(Y)))/nrow(X)
    Omega_1=(t(Y-X%*%B0)%*%(Y-X%*%B0) + lambda3*diag(ncol(Y)))/nrow(Y)
    invSigma1=solve(Sigma_1)
    invOmega1=solve(Omega_1)
    P1=t(chol(lambda1*Omega_1 +  lambda2*Sigma_1)) 
  
    #----------------------
    SVD_P1=svd(P1)
    U2_1=SVD_P1$u
    Sigma2_1=diag(SVD_P1$d)
    Sigma2vect_1=c(SVD_P1$d)
    tV2_1=SVD_P1$v
  
    #----------------------
    Btilde1= matrix(0, nrow=nrow(X), ncol=ncol(Y))
    S1=tV1%*%t(X)%*%Y%*%U2_1
    diagSigma1=(Sigma1vect)^2 #length p
    diagSigma2_1=(Sigma2vect_1)^2 #length s
    for(i in 1:(nrow(X))){
     for(j in 1:(ncol(Y))){
       Btilde1[i,j]<-S1[i,j]/(diagSigma1[i] + diagSigma2_1[j])
      }
    }
    B1=t(tV1)%*%Btilde1%*%t(U2_1)
  
    residual.cov_1=crossprod(Y-X%*%B1)/nrow(Y)
    obj_1=matrix.trace(residual.cov_1%*%invOmega1)-nrow(Y)*determinant(invOmega1, logarithm=TRUE)$mod[1] +
    lambda1*matrix.trace(t(B1)%*%B1) + lambda2*matrix.trace(invSigma1%*%t(B1)%*%B1) -ncol(X)*determinant(invSigma1, logarithm=TRUE)$mod[1] +
    lambda3*matrix.trace(invOmega1) +  lambda4*matrix.trace(invSigma1)
    #----------------------
    
    #----------------------
    tcontjmor <- sum(tcontjmor, 1)
    
    print(tcontjmor)
    
    
  }
  return(list(B1, invOmega1, invSigma1, tcontjmor))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in the MOR
CV_cmor<-function(lambda1, lambda2, lambda3, lambda4, X, Y, kfold=5) {
  lambda1 = sort(lambda1)
  lambda2 = sort(lambda2)
  lambda3 = sort(lambda3)
  lambda4 = sort(lambda4)
  n = nrow(Y)
  listt=list(lambda1, lambda2, lambda3,lambda4)
  dimll=length(lambda1)*length(lambda2)*length(lambda3)*length(lambda4)
  listindex=rep(list(rep(0,length(listt))),dimll )
  
  #CV_errors = array(0, c(dimll, kfold))
  CV_errors=rep(0,length(listindex))
  ind=sample(n)
  for (k in 1:kfold){
    
        leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    
    # validation set
    Y.valid = Y[leave.out,,drop=FALSE ]
    X.valid = X[leave.out,,drop=FALSE ]
    
    #this centers the training and validation data, uncomment if needed
    #Y.train_bar = apply(Y.train, 2, mean)
    #Y.train = scale(Y.train, center = Y.train_bar, scale = FALSE)
    #Y.valid = scale(Y.valid, center = Y.train_bar, scale = FALSE)
    
    #X.train_bar = apply(X.train, 2, mean)
    #X.train = scale(X.train, center = X.train_bar, scale = FALSE)
    #X.valid = scale(X.valid, center = X.train_bar, scale = FALSE)
    
    #X.train = scale(X.train, center = TRUE, scale = TRUE)
    #X.valid = scale(X.valid, attr(X.train, "scaled:center"), attr(X.train, "scaled:scale"))
    
    # loop over all tuning parameters
    
    for (m in 1:dimll) {
      for(i in 1:length(lambda1)) {
        for(j in 1:length(lambda2)) {
          for (t in 1:length(lambda3)) {
            for (l in 1:length(lambda4)) {
              # compute the joint penalized regression matrix estimate with the training set
              Estimes.train = cmor(lambda1[i], lambda2[j], lambda3[t], lambda4[l], X.train, Y.train) 
              B_hat.train = Estimes.train[[1]]
              invOmega_hat_train =  Estimes.train[[2]]
              invSigma_hat_train =  Estimes.train[[3]]
              listindex[[m]]=c(lambda1[i], lambda2[j], lambda3[t], lambda4[l])
              CV_errors[m]= CV_errors[m] + mean((Y.valid-X.valid%*%B_hat.train)^2)
             }
           }
         }
       }
     }   
   }
  
  # determine optimal tuning parameters
  AVG = mean(CV_errors) 
  bestindexAVG=which.min(CV_errors)
  error = min(CV_errors)
  opt = listindex[[bestindexAVG]]
  opt.lam1 = opt[1]
  opt.lam2 = opt[2]
  opt.lam3 = opt[3]
  opt.lam4 = opt[4]
  
  
  # return best lambdaO, lambdaB and other meaningful values
  return(list(lambda1=opt.lam1, lambda2=opt.lam2, lambda3=opt.lam3, lambda4=opt.lam4, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple ridge for B
multi.ridge <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=0,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple ridge
multi.ridge.cv <-function(X,Y, lam.vec, kfold=5, silent=TRUE){
  n=dim(Y)[1]
	q=dim(Y)[2]
  ind = sample(1:n) 
  p=dim(X)[2]
	err.vec=rep(0, length(lam.vec))
	
  for (k in 1:kfold)
  {
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
	  X.tr=X[-foldind, ]
	  Y.tr=Y[-foldind, ]
	  X.te=X[foldind, ]
	  Y.te=Y[foldind, ]
	  n.tr=dim(X.tr)[1]
	  
	  for(i in 1:length(lam.vec)){  
	    lam=lam.vec[i]
	  	bhatk = matrix(0, nrow=p, ncol=q)
	    for(kk in 1:q){
  	    bhatk[,kk]=as.numeric(glmnet(x=X.tr, y=Y.tr[,kk], family="gaussian", alpha=0,lambda=lam, standardize=FALSE)$beta)
		  }
		  err.vec[i]=err.vec[i]+mean((Y.te-X.te%*%bhatk)^2)  
    }
	  if(!silent) cat("Finished fold k = ", k, "\n")
	}
  best.B = matrix(0, nrow=p, ncol=q)
	best.i = which.min(err.vec)
	lam=lam.vec[best.i]
	for(kk in 1:q)
	{
	  best.B[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=0,lambda=lam, standardize=FALSE)$beta)
	}
  return(list(Bhat=best.B, lambda=lam.vec[best.i], cv.err=err.vec))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple lasso for B
multi.lasso <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=1,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple ridge
multi.lasso.cv <-function(X,Y, lam.vec, kfold=5, silent=TRUE){
  n=dim(Y)[1]
	q=dim(Y)[2]
  ind = sample(1:n) 
  p=dim(X)[2]
	err.vec=rep(0, length(lam.vec))
	
  for (k in 1:kfold)
  {
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
	  X.tr=X[-foldind, ]
	  Y.tr=Y[-foldind, ]
	  X.te=X[foldind, ]
	  Y.te=Y[foldind, ]
	
	  n.tr=dim(X.tr)[1]
	  for(i in 1:length(lam.vec))
    {  
	    lam=lam.vec[i]
	  	bhatk = matrix(0, nrow=p, ncol=q)
	    for(kk in 1:q)
  	  {
  	    bhatk[,kk]=as.numeric(glmnet(x=X.tr, y=Y.tr[,kk], family="gaussian", alpha=1,lambda=lam, standardize=FALSE)$beta)
		  }
		  err.vec[i]=err.vec[i]+mean((Y.te-X.te%*%bhatk)^2)  
    }
	  if(!silent) cat("Finished fold k = ", k, "\n")
	}
  best.B = matrix(0, nrow=p, ncol=q)
	best.i = which.min(err.vec)
	lam=lam.vec[best.i]
	for(kk in 1:q)
	{
	  best.B[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=1,lambda=lam, standardize=FALSE)$beta)
	}
  return(list(Bhat=best.B, lambda=lam.vec[best.i], cv.err=err.vec))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple ridge for B
multi.elasnet <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=.5,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple ridge
multi.elasnet.cv <-function(X,Y, lam.vec, kfold=5, silent=TRUE){
  n=dim(Y)[1]
	q=dim(Y)[2]
  ind = sample(1:n) 
  p=dim(X)[2]
	err.vec=rep(0, length(lam.vec))
	
  for (k in 1:kfold)
  {
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
	  X.tr=X[-foldind, ]
	  Y.tr=Y[-foldind, ]
	  X.te=X[foldind, ]
	  Y.te=Y[foldind, ]
	
	  n.tr=dim(X.tr)[1]
	  for(i in 1:length(lam.vec))
    {  
	    lam=lam.vec[i]
	  	bhatk = matrix(0, nrow=p, ncol=q)
	    for(kk in 1:q)
  	  {
  	    bhatk[,kk]=as.numeric(glmnet(x=X.tr, y=Y.tr[,kk], family="gaussian", alpha=.5,lambda=lam, standardize=FALSE)$beta)
		  }
		  err.vec[i]=err.vec[i]+mean((Y.te-X.te%*%bhatk)^2)  
    }
	  if(!silent) cat("Finished fold k = ", k, "\n")
	}
  best.B = matrix(0, nrow=p, ncol=q)
	best.i = which.min(err.vec)
	lam=lam.vec[best.i]

	for(kk in 1:q)
	{
	  best.B[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=.5,lambda=lam, standardize=FALSE)$beta)
	}
  return(list(Bhat=best.B, lambda=lam.vec[best.i], cv.err=err.vec))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the joint estimation for B and Omega 
JointEstim<-function(lambdaO, lambdaB, X, Y){

  # t=0 initialization
  diag_C_0=c(rep(1, nrow(t(X))))                    # Initialize C as the identity matrix of dim the number nrow of W_true
  invdiag_C_0=1/diag_C_0  
  B_ridge= ridgeEst(lambdaB, X, Y)
  Omega_ridge=t(t(Y)-t(B_ridge)%*%t(X))%*%(t(Y)-t(B_ridge)%*%t(X))+lambdaB*diag(ncol(t(Y)))
   
  Omega_0=solve(Omega_ridge)                               #initialize the inv Cov as omega ridge 
  B_0= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_0*t(X)%*%Omega_0%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_0*t(X)%*%Omega_0)))%*%X%*%(invdiag_C_0*t(X)%*%Omega_0)%*%Y), 8)
#----------------------
#update t=1 
    l21B_1=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_1[i]=2*norm((B_0[i,]),type = "2")
  }
  invdiag_C_1=l21B_1
  S_1=(1/nrow(t(Y)))*t(t(Y)-t(B_0)%*%t(X))%*%(t(Y)-t(B_0)%*%t(X))                                      
  Omega_1=glasso(S_1, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi         #get the precision matrix at t=1 from glasso
  B_1= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_1*t(X)%*%Omega_1%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_1*t(X)%*%Omega_1)))%*%X%*%(invdiag_C_1*t(X)%*%Omega_1)%*%Y), 8)
#-----------------------                   

#update t=2 
    l21B_2=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_2[i]=2*norm((B_1[i,]),type = "2")
  }
  invdiag_C_2=l21B_2
  S_2=(1/nrow(t(Y)))*t(t(Y)-t(B_1)%*%t(X))%*%(t(Y)-t(B_1)%*%t(X))                                      
  Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi         #get the precision matrix at t=1 from glasso
  B_2= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y), 8)
#----------------------- 
#start iteration in Algo 3
  tcontj=1

while((sum(l21B_2)<=sum(l21B_1))==TRUE  && tcontj<=itermax){  
#while(( (abs(sum(B_2 - B_1))>eps*abs(sum(B_ridge)))) && (tcontj<=itermax)){   #A different convergence criteria
  l21B_1=l21B_2
  invdiag_C_1=invdiag_C_2
  S_1=S_2
  Omega_1=Omega_2
  B_1=B_2                                             

#----------------------
  l21B_2=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_2[i]=2*norm((B_1[i,]),type = "2")
  }
  invdiag_C_2=l21B_2
  S_2=(1/nrow(t(Y)))*t(t(Y)-t(B_1)%*%t(X))%*%(t(Y)-t(B_1)%*%t(X))                                      
  Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi         #get the precision matrix at t=1 from glasso
  B_2= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y), 8)
  tcontj <- sum(tcontj, 1)                      
#---------------------- 
  print(tcontj)
                               
  }
  return(list(B_1, Omega_1, tcontj))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in the joint estimation Algorithm 3
CV_joint<-function(lambdaO, lambdaB, X, Y, kfold=5) {
  lambdaO = sort(lambdaO)
  lambdaB = sort(lambdaB)
  n = nrow(Y)
  CV_errors = matrix(0, nrow=length(lambdaO), ncol=length(lambdaB))
  ind=sample(n)
  for (k in 1:kfold){

    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    
    # validation set
    Y.valid = Y[leave.out,,drop=FALSE ]
    X.valid = X[leave.out,,drop=FALSE ]
    
    #this centers the training and validation data, uncomment if needed
    #Y.train_bar = apply(Y.train, 2, mean)
    #Y.train = scale(Y.train, center = Y.train_bar, scale = FALSE)
    #Y.valid = scale(Y.valid, center = Y.train_bar, scale = FALSE)
    
    #X.train_bar = apply(X.train, 2, mean)
    #X.train = scale(X.train, center = X.train_bar, scale = FALSE)
    #X.valid = scale(X.valid, center = X.train_bar, scale = FALSE)
    #X.train = scale(X.train, center = TRUE, scale = TRUE)
    #X.valid = scale(X.valid, attr(X.train, "scaled:center"), attr(X.train, "scaled:scale"))
              
    # loop over all tuning parameters
    for(i in 1:length(lambdaO)) {
	     for(j in 1:length(lambdaB)) {
	        # compute the joint penalized regression matrix estimate with the training set
                Estimes.train = JointEstim(lambdaO[i], lambdaB[j], X.train, Y.train) 
                B_hat.train = Estimes.train[[1]]
                Omega_hat_train =  Estimes.train[[2]]
                CV_errors[i,j] = CV_errors[i,j] + mean((Y.valid-X.valid%*%B_hat.train)^2)
             }
         }
    }
    
    # determine optimal tuning parameters
    AVG = apply(CV_errors, 1, mean)
    error = min(AVG)
    opt = which.min(CV_errors) %% (dim(CV_errors)[1])
    opt = (opt != 0)*opt + (opt == 0)*(dim(CV_errors)[1])
    opt.i = opt
    opt.j = which.min(CV_errors[opt,])
    opt.lamO = lambdaO[opt.i]
    opt.lamB = lambdaB[opt.j]
  
  
    # return best lambdaO, lambdaB and other meaningful values
  return(list(lambdaO=opt.lamO, lambdaB=opt.lamB, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Synthetic data generation 
itermax=100
n=70              			# nrow_Y=n and sample size
p=800   		        	# nrow_B/ncol_X
s=30   		        		# ncol_B/Y
m=p+n
rho_E=.9 		     		# normaly =c(0,.4,.9)
s_1=.5   		    		# Probability of success for independent entries in K
s_2=.4   		    		# Probability of success for the ones in Q
reps=20					    # number of replicates
eps=10^-5
n.val=50

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# '@ generate data for simulations
# 'replicate of predictor matrices X
listX=list()
for (l in 1:reps){
  SigmaX = matrix(0.7, nrow = p, ncol = p)
  for (i in 1:p){
    for (j in 1:p){
      SigmaX[i, j] = SigmaX[i, j]^abs(i - j)
    }
  }
  Z = matrix(rnorm(n*p), nrow = n, ncol = p)
  out = eigen(SigmaX, symmetric = TRUE)
  SigmaX.sqrt = out$vectors %*% diag(out$values^0.5)
  SigmaX.sqrt = SigmaX.sqrt %*% t(out$vectors)
  X = Z %*% SigmaX.sqrt
  listX[[l]]= X  
}

#-----------------------------------------------------------------------------------
# 'replicate of B matrices 

listB=list()
for (l in 1:reps){
  W=matrix(rnorm(p*s,0,1),p,s)
  K=matrix(rbinom(p*s,1,s_1),p,s)
  Q=matrix(0,p,s)
  q=rbinom(p,1,s_2)
  for (i in 1:p){
    if (q[i]>0)
      Q[i,]=rep(1,s)
    
  }
  #B=W*K*Q
  #the other way we can simulate a sparse matrix B is as follows
  #listB[[l]]=matrix(rbinom(p*s, size=1, prob=s_1)*runif(p*s, min=1, max=2), nrow=p, ncol=s)
  listB[[l]]=W*Q + K*W
}

#-----------------------------------------------------------------------------------
# 'replicate of error matrices 
listE=list()
for (l in 1:reps){
  SigmaE = matrix(rho_E, nrow = s, ncol = s)
  for (i in 1:s){
    for (j in 1:s){
      SigmaE[i, j] = SigmaE[i, j]^abs(i - j)
    }
  }
  set.seed(20)
  P = matrix(rnorm(n*s), nrow = n, ncol = s)
  out = eigen(SigmaE, symmetric = TRUE)
  SigmaE.sqrt = out$vectors %*% diag(out$values^0.5)
  SigmaE.sqrt = SigmaE.sqrt %*% t(out$vectors)
  E = P %*% SigmaE.sqrt
  listE[[l]]= E   
}

#-----------------------------------------------------------------------------------
# 'replicate of Y matrices 
listY=list()
for (l in 1:reps){
  listY[[l]]=listX[[l]]%*%listB[[l]] + listE[[l]]
}
#-----------------------------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 1 l21 features selection with CV
# The l21fs estimate with CV
listX.val.l21fs=list()
listY.val.l21fs=list()
listX.tr.l21fs=list()
listY.tr.l21fs=list()
Bhat_l21fs_Sim = list()
lambdaB.opt.l21fs=c()
for (l in 1:reps) {

  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.l21fs[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.l21fs[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.l21fs[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.l21fs[[l]]=listY[[l]][-c(1:n.val),]
  
  # 'find the optimal lambdaB with 5-folds CV
  L.opt=CV_feat_sel(listY.tr.l21fs[[l]], listX.tr.l21fs[[l]], lam=seq(8,12,.5), kfold=5)
  lambdaB.opt.l21fs[l]=L.opt[[1]]
  Bhat_l21fs_Sim[[l]] = FeatSelect(lambdaB.opt.l21fs[l], listX.tr.l21fs[[l]], listY.tr.l21fs[[l]])
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 2 ridge Estimate with CV using glmnet package
listX.val.multiridge=list()
listY.val.multiridge=list()
listX.tr.multiridge=list()
listY.tr.multiridge=list()
listY.tr.multiridge_bar = list()
Bhat_multiridge_Sim = list()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.multiridge[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.multiridge[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.multiridge[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.multiridge[[l]]=listY[[l]][-c(1:n.val),]
  
  # 'find the optimal lambdaB with 5-folds CV
  L.opt.multiridge=multi.ridge.cv(listX.tr.multiridge[[l]],listY.tr.multiridge[[l]], lam.vec=2^seq(-3, 3, 1), kfold=5, silent=TRUE)

  Bhat_multiridge_Sim[[l]] = L.opt.multiridge[[1]]
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 3 multi Lasso Estimate with CV using glmnet package
listX.val.multilasso=list()
listY.val.multilasso=list()
listX.tr.multilasso=list()
listY.tr.multilasso=list()
Bhat_multilasso_Sim = list()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.multilasso[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.multilasso[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.multilasso[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.multilasso[[l]]=listY[[l]][-c(1:n.val),]
 
  # 'find the optimal lambdaB with 5-folds CV
  L.opt.multilasso=multi.lasso.cv(listX.tr.multilasso[[l]],listY.tr.multilasso[[l]], lam.vec=2^seq(-3, 3, 1), kfold=5, silent=TRUE)

  Bhat_multilasso_Sim[[l]] = L.opt.multilasso[[1]]
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 4 Elastic net Estimate with CV using glmnet package
listX.val.multielasnet=list()
listY.val.multielasnet=list()
listX.tr.multielasnet=list()
listY.tr.multielasnet=list()
Bhat_multielasnet_Sim = list()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.multielasnet[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.multielasnet[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.multielasnet[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.multielasnet[[l]]=listY[[l]][-c(1:n.val),]
  
  # 'find the optimal lambdaB with 5-folds CV
  L.opt.multielasnet=multi.elasnet.cv(listX.tr.multielasnet[[l]],listY.tr.multielasnet[[l]], lam.vec=2^seq(-3, 3, 1), kfold=5, silent=TRUE)

  Bhat_multielasnet_Sim[[l]] = L.opt.multielasnet[[1]]
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 5 cMOR with CV
#joint estim with CV on (lamB, lamO) and covariance estim
listX.val.cmor=list()
listY.val.cmor=list()
listX.tr.cmor=list()
listY.tr.cmor=list()
Bhat_cmor_Sim=list()
invOmega_cmor_Sim=list()
invSigma_cmor=list()
lam1.opt.cmor=c()
lam2.opt.cmor=c()
lam3.opt.cmor=c()
lam4.opt.cmor=c()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.cmor[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.cmor[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.cmor[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.cmor[[l]]=listY[[l]][-c(1:n.val),]

  L.cmor.opt=CV_cmor(lambda1 = 2^seq(-2,-1,1), lambda2 = 2^seq(-2,-1,1), lambda3 = 2^seq(-2,-1,1), lambda4 = 2^seq(-2,-1,1), listX.tr.cmor[[l]], listY.tr.cmor[[l]], kfold=5)
  lam1.opt.cmor[l]=L.cmor.opt[[1]]
  lam2.opt.cmor[l]=L.cmor.opt[[2]]
  lam3.opt.cmor[l]=L.cmor.opt[[3]]
  lam4.opt.cmor[l]=L.cmor.opt[[4]]

  Final_cmor_est_CV=cmor(lam1.opt.cmor[l], lam2.opt.cmor[l], lam3.opt.cmor[l], lam4.opt.cmor[l], listX.tr.cmor[[l]], listY.tr.cmor[[l]])
  Bhat_cmor_Sim[[l]]=Final_cmor_est_CV[[1]]
  invOmega_cmor_Sim[[l]]=Final_cmor_est_CV[[2]]
  invSigma_cmor[[l]]=Final_cmor_est_CV[[3]]

}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 6 remMap with features selection and CV
#joint estim with CV on (lamB, lamO) and covariance estim
listX.val.remMap=list()
listY.val.remMap=list()
listX.tr.remMap=list()
listY.tr.remMap=list()
Bhat_remMap_Sim=list()
lamL1.pick=c()
lamL2.pick=c()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.remMap[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.remMap[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.remMap[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.remMap[[l]]=listY[[l]][-c(1:n.val),]

lamL1.v=2^seq(-2,-1,1)#2^seq(-5,-2,1)
lamL2.v=2^seq(-2,-1,1)#2^seq(-5,-2,1)
L.remMap=remMap.CV(X=listX.tr.remMap[[l]], Y=listY.tr.remMap[[l]], lamL1.v, lamL2.v, C.m=NULL, fold=5, seed=123)
pick=which.min(as.vector(L.remMap$ols.cv))
lamL1.pick[l]=L.remMap$l.index[1,pick]    ##find the optimal (LamL1,LamL2) based on the cv score
lamL2.pick[l]=L.remMap$l.index[2,pick]
##fit the remMap model under the optimal (LamL1,LamL2).
Final_remMAP=remMap(listX.tr.remMap[[l]], listY.tr.remMap[[l]],lamL1=lamL1.pick[l], lamL2=lamL2.pick[l], phi0=NULL, C.m=NULL)
print(paste("lamL1=", round(lamL1.pick[l],3), "; lamL2=", round(lamL2.pick[l],3), sep=""))
Bhat_remMap_Sim[[l]]=Final_remMAP$phi

}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 7 Multitrait Bayes B
#joint estim with CV on (lamB, lamO) and covariance estim
listX.val.MBayesB=list()
listY.val.MBayesB=list()
listX.tr.MBayesB=list()
listY.tr.MBayesB=list()
SVD=list()
U=list()
D=list()
V=list()
BETA_MBayesB_Sim=list()
Bhat_MBayesB_Sim=list()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.MBayesB[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.MBayesB[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.MBayesB[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.MBayesB[[l]]=listY[[l]][-c(1:n.val),]
  
  # Singular value decompositions Y=UDV'
  SVD[[l]]=svd(listY.tr.MBayesB[[l]])
  U[[l]]=SVD[[l]]$u
  D[[l]]=diag(SVD[[l]]$d)
  V[[l]]=SVD[[l]]$v
  
  
  Bhat_MBayesB_Sim[[l]]=matrix(nrow=ncol(listX.tr.MBayesB[[l]]),ncol=ncol(listY.tr.MBayesB[[l]]))
  
  ETA=list(list(X=listX.tr.MBayesB[[l]],model='BayesB'))
  for(i in 1:ncol(listY.tr.MBayesB[[l]])){
    fm=BGLR(y=U[[l]][,i],ETA=ETA,verbose=F) #use more iterations! if not converge. Here we can also set the Pi values if needed
    Bhat_MBayesB_Sim[[l]][,i]=fm$ETA[[1]]$b
  }
  
  # Rotating coefficients to put them in marker space
  BETA_MBayesB_Sim[[l]]=Bhat_MBayesB_Sim[[l]]%*%D[[l]]%*%t(SVD[[l]]$v)
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 8 jointl21 with features selection and CV
#joint estim with CV on (lamB, lamO) and covariance estim
listX.val.l21joint=list()
listY.val.l21joint=list()
listX.tr.l21joint=list()
listY.tr.l21joint=list()

Bhat_l21joint_Sim=list()
Omega_hat_l21joint_Sim=list()
lamO.opt.l21joint=c()
lamB.opt.l21joint=c()
for (l in 1:reps) {
  
  listY[[l]]=scale(listY[[l]], scale = TRUE) 
  listX[[l]]=scale(listX[[l]], scale = TRUE)
  
  listY.tr.l21joint[[l]]=listY[[l]][c(1:n.val),]
  listX.tr.l21joint[[l]]=listX[[l]][c(1:n.val),]
  
  listX.val.l21joint[[l]]=listX[[l]][-c(1:n.val),]
  listY.val.l21joint[[l]]=listY[[l]][-c(1:n.val),]
  
  L.l21joint=CV_joint(lambdaO=2^seq(-11, -9, 1), lambdaB=seq(.3, .5, .1), listX.tr.l21joint[[l]], listY.tr.l21joint[[l]], kfold=5)
  lamO.opt.l21joint[l]=L.l21joint[[1]]
  lamB.opt.l21joint[l]=L.l21joint[[2]]
  
  Final_joint_est_CVfs=JointEstim(lamO.opt.l21joint[l], lamB.opt.l21joint[l], listX.tr.l21joint[[l]], listY.tr.l21joint[[l]]) 
  
  Bhat_l21joint_Sim[[l]] = Final_joint_est_CVfs[[1]]
  Omega_hat_l21joint_Sim[[l]]=Final_joint_est_CVfs[[2]]
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
