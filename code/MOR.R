# Function to implement GS using multiple output regression
# Y nxs response matrix
# X nxp predictor matrix 
# lambda1, lambda2, lambda3 and lambda4 are the regularization parameters, see (He et all. 2016) for details.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MOR <- function(lambda1, lambda2, lambda3, lambda4, X, Y){
  #----------------------
  
  # t=0 initialization
  
  Sigma_0=diag(ncol(Y))                           # Initialize Omega as the identity matrix of dim the number ncol of Y 
  Omega_0=diag(ncol(Y))                           # Initialize Sigma as the identity matrix of dim the number ncol of Y 
  invSigma0=solve(Sigma_0)
  invOmega0=solve(Omega_0)
  P0=t(chol(lambda1*Omega_0 +  lambda2*Sigma_0))  # to get the lower triangular matrix P (notice the transpose, otherwise the upper triangular matrix is provided).
  
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
  V2_1=SVD_P1$v
  
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
  
  #----------------------
  
  tcontjmor=1
  while((abs(sum(B1))<=abs(sum(B0)))==TRUE  && tcontjmor<=itermax){  
    
    Sigma_0=Sigma_1
    Omega_0=Omega_1  
    invSigma0=invSigma1                                         
    invOmega0=invOmega1   
    P0=P1
    #----------------------
    SVD_P0=svd(P0)
    U2_0=SVD_P0$u
    Sigma2_0=diag(SVD_P0$d)
    Sigma2vect_0=c(SVD_P0$d)
    V2_0=SVD_P0$v
    
    #---------------------- 
    Btilde0= matrix(0, nrow=nrow(X), ncol=ncol(Y))
    S0=tV1%*%t(X)%*%Y%*%U2_0
    diagSigma1=(Sigma1vect)^2          #length p
    diagSigma2_0=(Sigma2vect_0)^2      #length s
    for(i in 1:(nrow(X))){
      for(j in 1:(ncol(Y))){
        Btilde0[i,j]<-S0[i,j]/(diagSigma1[i] + diagSigma2_0[j])
      }
    }
    B0=t(tV1)%*%Btilde0%*%t(U2_0)
    
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
    V2_1=SVD_P1$v
    
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
    
    #----------------------
    tcontjmor <- sum(tcontjmor, 1)
    
    print(tcontjmor)
    
    
  }
  return(list(B1, invOmega1, invSigma1))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
