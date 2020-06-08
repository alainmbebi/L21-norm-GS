# Function to implement Cross Validation for L21_joint_estim
# Y nxs response matrix
# X nxp predictor matrix 
# lambdaO vector of positive tuning parameters for the precision matrix, default grid is 2^seq(-3, 3, 1).
# lambdaB vector of positive tuning parameters for the regression coefficient matrix, default grid is 2^seq(-3, 3, 1).
# kfold is the number of folds for cross validation.
# The returned list includes
# lambdaO--> optimal tuning parameter for the precision matrix.
# lambdaB--> optimal tuning parameter for the regression coefficient matrix.
# min.error--> minimum average cross validation error for optimal parameters.
# avg.error--> average cross validation error across all folds.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CV_L21_joint_estim <- function(lambdaO, lambdaB, X, Y, kfold=5) {
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
    Y.valid = Y[leave.out, ]
    X.valid = X[leave.out, ]
    
    # loop over all tuning parameters
    for(i in 1:length(lambdaO)) {
      for(j in 1:length(lambdaB)) {
        # compute the joint penalized regression matrix estimate with the training set
        Estimes.train = L21_joint_estim(lambdaO[i], lambdaB[j], X.train, Y.train) 
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
