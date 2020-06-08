#function for the joint estimation for the regression coefficient (B) and precision (Omega) matrices 
# Y nxs response matrix
# X nxp predictor matrix 
# lambdaB regularization parameter for the regression coefficient matrix B
# lambdaO regularization parameter for the precision matrix \Omega
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

L21_joint_estim <- function(lambdaO, lambdaB, X, Y){
  #----------------------
  # t=0 initialization
  diag_C_0=c(rep(1, nrow(t(X))))            
  invdiag_C_0=1/diag_C_0  
  B_Ridge= Ridge_estim(lambdaB, X, Y)
  Omega_ridge=t(t(Y)-t(B_Ridge)%*%t(X))%*%(t(Y)-t(B_Ridge)%*%t(X))+lambdaB*diag(ncol(t(Y)))
  
  Omega_0=solve(Omega_ridge)                               #initialize the inv Cov as omega ridge 
  B_0= (2/(nrow(t(Y))*lambdaB))*invdiag_C_0*t(X)%*%Omega_0%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_0*t(X)%*%Omega_0)))%*%X%*%(invdiag_C_0*t(X)%*%Omega_0)%*%Y)
  #----------------------
 
  #update t=1 
  l21B_1=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_1[i]=2*norm((B_0[i,]),type = "2")
  }
  invdiag_C_1=l21B_1
  S_1=(1/nrow(t(Y)))*t(t(Y)-t(B_0)%*%t(X))%*%(t(Y)-t(B_0)%*%t(X))                                      
  Omega_1=glasso(S_1, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi         #get the precision matrix from glasso
  B_1= (2/(nrow(t(Y))*lambdaB))*invdiag_C_1*t(X)%*%Omega_1%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_1*t(X)%*%Omega_1)))%*%X%*%(invdiag_C_1*t(X)%*%Omega_1)%*%Y)
  #-----------------------                   
  
  #update t=2 
  l21B_2=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_2[i]=2*norm((B_1[i,]),type = "2")
  }
  invdiag_C_2=l21B_2
  S_2=(1/nrow(t(Y)))*t(t(Y)-t(B_1)%*%t(X))%*%(t(Y)-t(B_1)%*%t(X))                                      
  Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi      
  B_2= (2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y)
  #----------------------- 
  #start iteration
 
  tcontj=1
  #while((sum(l21B_2)<=sum(l21B_1))==TRUE  && tcontj<=itermax){ ## other convergence criteria
  while(( (abs(sum(B_2 - B_1))>eps*abs(sum(B_Ridge)))) && (tcontj<=itermax)){ 
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
    Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi      
    B_2= (2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y)
    tcontj <- sum(tcontj, 1)                      
    #---------------------- 
    print(tcontj)
    
  }
  return(list(B_1, Omega_1))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
