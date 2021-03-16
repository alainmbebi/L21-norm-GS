
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #start
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set.seed(143)

#load required libraries
library(MASS)
library(MRCE)       # to implement MRCE
library(glasso)     # to implement glasso
library(glmnet)     # to implement the multivariate LASSO
library(MTGS)       # to implement all MTGS models and also provides the Brassica napus data 
library(BGLR)       # to implement GBLUP, MBayesB and also provides the wheat data 
library(FactoMineR) # to compute the RV coefficient and test its significance
library(hydroGOF)   # contains various function for summary statistics, such as mse, mae...

#source the required functions
source("CV_L21_featselect.R")
source("CV_L21_joint_estim.R")
source("CV_MOR.R")
source("CV_Ridge.R")
source("L21_featselect.R")
source("L21_joint_estim.R")
source("MOR.R")
source("Ridge_estim.R")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
itermax=200  		    # maximum iterations
eps=10^-5             #threshold for convergence in L21-joint as described in the manuscript
data(brassica_data)
X<-as.matrix(brassica_data[,1:100])
Y<-as.matrix(brassica_data[,101:103])
Y.tr=Y[c(1:40),]    #training population (phenotypes)
X.tr=X[(1:40),]     #training population (genotypes)
Y.hold=Y[-c(1:40),] #selection candidate (phenotypes)
X.hold=X[-c(1:40),] #selection candidates (genotypes)
X.tr = scale(X.tr, center = TRUE, scale = TRUE)
X.hold <-scale(X.hold, attr(X.tr, "scaled:center"), attr(X.tr, "scaled:scale"))
r<-0.35 # see the MTGS package for detail

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 1: Ridge Estimate
# the Ridge estimate
# find the optimal lambdaB with 5-folds CV
L.opt.Ridge=CV_Ridge(Y.tr, X.tr, lambdaB=2^seq(-2,2,1), kfold=5)
lambdaB.opt.Ridge=L.opt.Ridge[[1]]    #the optimal value was 4

Bhat_Ridge= Ridge_estim(lambdaB.opt.Ridge, X.tr, Y.tr)
Pred_Ridge=X.hold%*%Bhat_Ridge

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 2: l21 features selection with CV
# find the optimal lambdaB with 5-folds CV
L.opt=CV_L21_featselect(Y.tr, X.tr, lam=2^seq(-4,-2,1), kfold=5)
lambdaB.opt.l21fs=L.opt[[1]]    #the optimal value was .25 

Bhat_l21fs = L21_featselect(lambdaB.opt.l21fs, X.tr, Y.tr)
Pred_l21fs=X.hold%*%Bhat_l21fs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 3: cMOR
# find the optimal lambdas with 5-folds CV
L.MOR.opt=CV_MOR(lambda1 = 2^seq(-2,-1,1), lambda2 = 2^seq(-2,-1,1), lambda3 = 2^seq(-2,-1,1), lambda4 = 2^seq(-2,-1,1), X.tr, Y.tr,  kfold=5)
lam1.opt.MOR=L.MOR.opt[[1]]   #the optimal value was .5 
lam2.opt.MOR=L.MOR.opt[[2]]   #the optimal value was .5 
lam3.opt.MOR=L.MOR.opt[[3]]   #the optimal value was .5 
lam4.opt.MOR=L.MOR.opt[[4]]   #the optimal value was .5 

Final_MOR_est_CV=MOR(lam1.opt.MOR, lam2.opt.MOR, lam3.opt.MOR, lam4.opt.MOR, X.tr, Y.tr)
Bhat_MOR=round(Final_MOR_est_CV[[1]], 6)
invOmega_MOR=Final_MOR_est_CV[[2]]
invSigma_MOR=Final_MOR_est_CV[[3]]
Pred_MOR=X.hold%*%Bhat_MOR

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 4: mlasso
Pred_MTGS.mlasso.tr=MTGS.mlasso(X.tr,Y.tr,r)$Pred
Bhat_MTGS.mlasso=ginv(X.tr)%*%Pred_MTGS.mlasso.tr
Pred_MTGS.mlasso=X.hold%*%Bhat_MTGS.mlasso

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 5: MRCE
Bhat_MTGS.mrce=MTGS.mrce(X.tr,Y.tr,r)$Bhat
Pred_MTGS.mrce=X.hold%*%Bhat_MTGS.mrce

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 6: kmLASSO
Pred_MTGS.kmlasso.tr=MTGS.kmlasso(X.tr, Y.tr)$Pred[,,1]
Bhat_MTGS.kmlasso=ginv(X.tr)%*%Pred_MTGS.kmlasso.tr
Pred_MTGS.kmlasso=X.hold%*%Bhat_MTGS.kmlasso

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 7: GBLUP

ETA=list(list(X.tr=X.tr,model="BRR"))
fmR1<-BGLR(y=Y.tr[,1],ETA=ETA,nIter=1000,burnIn=500,thin=10)
fmR2<-BGLR(y=Y.tr[,2],ETA=ETA,nIter=1000,burnIn=500,thin=10)
fmR3<-BGLR(y=Y.tr[,3],ETA=ETA,nIter=1000,burnIn=500,thin=10)
bHat1_GBLUP<- fmR1$ETA[[1]]$b
bHat2_GBLUP<- fmR2$ETA[[1]]$b
bHat3_GBLUP<- fmR3$ETA[[1]]$b      
Bhat_GBLUP=cbind(bHat1_GBLUP,bHat2_GBLUP,bHat3_GBLUP)
Pred_GBLUP=X.hold%*%Bhat_GBLUP

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 8 Multitrait Bayes B
  SVD=svd(Y.tr)
  U=SVD$u
  D=diag(SVD$d)
  V=SVD$v
  
  
 Bhat_MBayesB=matrix(nrow=ncol(X.tr),ncol=ncol(Y.tr))

 ETA=list(list(X=X.tr,model='BayesB'))
 for(i in 1:ncol(Y.tr)){
	fm=BGLR(y=U[,i],ETA=ETA,verbose=F) #use more iterations if needed!
	Bhat_MBayesB[,i]=fm$ETA[[1]]$b
 }

 # Rotating coefficients to put them in marker space
BETA_MBayesB=Bhat_MBayesB%*%D%*%t(SVD$v)
Pred_MBayesB=X.hold%*%BETA_MBayesB

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#L21-joint estimation with CV 
L.joint.optfs=CV_L21_joint_estim(lambdaO=12, lambdaB=2^seq(-4,-3,1), X.tr, Y.tr,  kfold=5)
lamO.opt.jointfs=L.joint.optfs[[1]]   #the optimal value was 12 
lamB.opt.jointfs=L.joint.optfs[[2]]   #the optimal value was .125 

Final_joint_est_CVfs=L21_joint_estim(lamO.opt.jointfs, lamB.opt.jointfs, X.tr, Y.tr)
Bhat_rmrce_jointfs = round(Final_joint_est_CVfs[[1]],6)
Omega_hat_rmrce_jointfs=Final_joint_est_CVfs[[2]]
Pred_rmrce_jointfs=X.hold%*%Bhat_rmrce_jointfs


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                   #end
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
