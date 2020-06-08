#function for cross validation in the MOR

CV_MOR <- function(lambda1, lambda2, lambda3, lambda4, X, Y, kfold=5) {
  lambda1 = sort(lambda1)
  lambda2 = sort(lambda2)
  lambda3 = sort(lambda3)
  lambda4 = sort(lambda4)
  n = nrow(Y)
  listt=list(lambda1, lambda2, lambda3,lambda4)
  dimll=length(lambda1)*length(lambda2)*length(lambda3)*length(lambda4)
  listindex=rep(list(rep(0,length(listt))),dimll )
  
  CV_errors=rep(0,length(listindex))
  ind=sample(n)
  for (k in 1:kfold){
    
    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    #this centers the data uncomment if needed
    #Y_bar = apply(Y.train, 2, mean)
    #X_bar = apply(X.train, 2, mean)
    #Y.train = scale(Y.train, center = Y_bar, scale = FALSE)
    #X.train = scale(X.train, center = X_bar, scale = FALSE)
    
    # validation set
    Y.valid = Y[leave.out, ]
    X.valid = X[leave.out, ]
    #this centers the data uncomment if needed
    #Y.valid = scale(Y.valid, center = Y_bar, scale = FALSE)
    #X.valid = scale(X.valid, center = X_bar, scale = FALSE)
    
    # loop over all tuning parameters
    
    for (m in 1:dimll) {
      for(i in 1:length(lambda1)) {
        for(j in 1:length(lambda2)) {
          for (t in 1:length(lambda3)) {
            for (l in 1:length(lambda4)) {
              # compute the joint penalized regression matrix estimate with the training set
              Estimes.train = MOR(lambda1[i], lambda2[j], lambda3[t], lambda4[l], X.train, Y.train) 
              B_hat.train = Estimes.train[[1]]
              invOmega_hat_train =  Estimes.train[[2]]
              invSigma_hat_train =  Estimes.train[[3]]
              #for (s in 1:dimll) {
              listindex[[m]]=c(lambda1[i], lambda2[j], lambda3[t], lambda4[l])
              #CV_errors[m,k]= CV_errors[m,k] + mean((Y.valid-X.valid%*%B_hat.train)^2)
              CV_errors[m]= CV_errors[m] + mean((Y.valid-X.valid%*%B_hat.train)^2)
            }
          }
          
          #}
          
        }
      }
      #CV_errors[m]= CV_errors[m] + mean((Y.valid-X.valid%*%B_hat.train)^2)
    }
    
  }
  
  # determine optimal tuning parameters
  AVG = mean(CV_errors) 
  #AVG = mean(CV_errors)
  bestindexAVG=which.min(CV_errors)
  #error = lapply(AVG, min) 
  error = min(CV_errors)
  opt = listindex[[bestindexAVG]]
  opt.lam1 = opt[1]
  opt.lam2 = opt[2]
  opt.lam3 = opt[3]
  opt.lam4 = opt[4]
  
  
  # return best lambdaO, lambdaB and other meaningful values
  return(list(lambda1=opt.lam1, lambda2=opt.lam2, lambda3=opt.lam3, lambda4=opt.lam4, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}
