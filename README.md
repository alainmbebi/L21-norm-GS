# L21-norm-GS

R functions to implement genotype-phenotype association and genomic prediction using ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression, as described in [(A. Mbebi et all. 2020)]().

The following R scripts are contained in the 'code' folder:

* L21_featselect.R runs the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression as in [(Nie et all. 2010)](http://papers.nips.cc/paper/3988-efficient-and-robust-feature-selection-via-joint-l21-norms-minimization)

* CV_L21_featselect.R selects the tuning parameter for L21_featselect.R using K-folds cross-validation

* RidgeEst.R compute the Ridge estimate

* CV_Ridge.R selects the tuning parameter for RidgeEst.R using K-folds cross-validation

* L21_joint_estim.R performs the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression and jointly estimate the regression coefficient and precision matrices

* CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

* MOR.R runs multiple output regression [(He et all. 2016)](https://academic.oup.com/bioinformatics/article/32/12/i37/2288681)

* CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

* MTGS.mlasso.R is called from MTGS package and compute the GEBVs using multivariate LASSO 

* MTGS.kmlasso.R is called from MTGS package and compute the GEBVs using Kernelized multivariate LASSO

* MTGS.mrce.R is called from MTGS package and compute the GEBVs using MRCE

* coeffRV.R


By default the all cross-validation scripts use 5-folds

MTGS: Multi Traits Genomic Selection
GS: Genomic Selection
CV: Cross Validation
LASSO: Least Absolute Shrinkage and Selection Operator
MOR: Multiple Output Regression
MRCE: Multivariate Regression with Covariance Estimation

emma.r use the emma function (Kang et all. 2008) for the MRLE
