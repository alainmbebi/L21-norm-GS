# Function to implement Cross Validation for L21_featselect
# Y nxs response matrix
# X nxp predictor matrix 
# lam vector of positive tuning parameters, default grid is 2^seq(-3, 3, 1).
# kfold is the number of folds for cross validation.
# The returned list includes
# lam--> optimal tuning parameter.
# min.error--> minimum average cross validation error for optimal parameters.
# avg.error--> average cross validation error across all folds.


CV_L21_featselect <- function(Y, X, lam, kfold=5) {
  
  lam = sort(lam)
  # initialize
  CV_errors = array(0, c(length(lam), kfold))
  
  n = nrow(Y)
  ind = sample(n)
  
  # parsing data into folds and performing CV
  
  for (k in 1:kfold) {
    
    leave.out = ind[ (1 + floor((k - 1)*n/kfold)):floor(k*n/kfold) ]
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    
    # validation set
    Y.valid = Y[leave.out, ]
    X.valid = X[leave.out, ]
        
    # loop over all tuning parameters
    for (i in 1:length(lam)) {
      
      # compute the penalized regression matrix estimate with the training set
      B_hat.train=L21_featselect(lam[i], X.train, Y.train) 
      
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
