# Function to implement GS for L21-norm regularized multivariate regression using the model in (Nie et all. 2010)
# Y nxs response matrix
# X nxp predictor matrix 
# lambdaB regularization parameter
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

L21_featselect<-function(lambdaB, X, Y){
  #--------------------------
  A=cbind(X, lambdaB*diag(nrow(Y)))
  B_l21fs=matrix(0, nrow=ncol(X), ncol=ncol(Y))
  
  #--------------------------
  # initialization t=0 
  diag_D_0=c(rep(1, (ncol(X)+nrow(Y))))               
  invdiag_D_0=1/diag_D_0              			       
  
  #--------------------------
  #update of W at t=2
  
  W_1=round((invdiag_D_0*t(A))%*%solve(A%*%(invdiag_D_0*t(A)))%*%Y ,6)   
  l21W_1=c(rep(0, (ncol(X)+nrow(Y))))
  for(i in 1:(ncol(X)+nrow(Y))){
    l21W_1[i]=2*norm((W_1[i,]),type = "2")
  }

  invdiag_D_1=l21W_1
  
  #--------------------------
  #update of W at t=2
  
  W_2=round((invdiag_D_1*t(A))%*%solve(A%*%(invdiag_D_1*t(A)))%*%Y ,6) 
  l21W_2=c(rep(0, (ncol(X)+nrow(Y))))
  for(i in 1:(ncol(X)+nrow(Y))){
    l21W_2[i]=2*norm((W_2[i,]),type = "2")
  }

  invdiag_D_2=l21W_2
  
  #--------------------------
  
  tcont=0
  while((sum(l21W_2)<=sum(l21W_1))==TRUE  && tcont<=itermax){
    W_1=W_2
    l21W_1=l21W_2                                                                   
    invdiag_D_1=invdiag_D_2                                                         
    
    W_2=round((invdiag_D_1*t(A))%*%solve(A%*%(invdiag_D_1*t(A)))%*%Y ,6)        
    l21W_2=c(rep(0,(ncol(X)+nrow(Y))))
    for(i in 1:(ncol(X)+nrow(Y))){
      l21W_2[i]=2*norm((W_2[i,]),type = "2")
    }
    invdiag_D_2=l21W_2
    tcont <- sum(tcont, 1)
    
    print(tcont)
    B_l21fs=round(W_1[1:ncol(X),],6)      #extract the estimated regression coefficient matrix
    
  }
  return(B_l21fs)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
