#function for cross validation in elastic net 
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
