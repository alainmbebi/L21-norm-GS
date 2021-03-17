# Function to compute Ridge estimate
# Y nxs response matrix
# X nxp predictor matrix 
# lambdaB regularization parameter
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
