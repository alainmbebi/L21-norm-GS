# L21-norm-GS

R functions to implement genotype-phenotype association and genomic prediction using ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression, as described in [(A. Mbebi et all. 2020)]().

The following R scripts are contained in the 'code' folder:

* L21_featselect.R which uses the model proposed in [(Nie et all. 2010)](http://papers.nips.cc/paper/3988-efficient-and-robust-feature-selection-via-joint-l21-norms-minimization) to implement GS under the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression 

* CV_L21_featselect.R which selects the tuning parameter for L21_featselect.R using K-folds cross-validation

* Ridge_estim.R compute the Ridge estimate

* CV_Ridge.R selects the tuning parameter for Ridge_estim.R using K-folds cross-validation

* L21_joint_estim.R performs GS using the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression that jointly estimates the regression coefficient and precision matrices

* CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

* MOR.R runs multiple output regression [(He et all. 2016)](https://academic.oup.com/bioinformatics/article/32/12/i37/2288681)

* CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

* MTGS.mlasso.R uses the MTGS.mlasso function from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using multivariate LASSO 

* MTGS.kmlasso.R uses the MTGS.kmlasso function from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using Kernelized multivariate LASSO

* MTGS.mrce.R uses the MTGS.mrce function from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using MRCE

* coeffRV.R uses the coeffRV function from [FactoMineR R package](https://CRAN.R-project.org/package=FactoMineR) and is used to compute the RV coefficient

Depends: R (>= 3.6)
Imports: glmnet, kernlab, MRCE, BGLR, FactoMineR
Licence: GPL-3
By default, all cross-validation scripts use 5-folds

Abreviations

CV: Cross Validation, GBLUP: Genomic Best Linear Unbiased Prediction, GS: Genomic Selection, 
LASSO: Least Absolute Shrinkage and Selection Operator, MOR: Multiple Output Regression
MRCE: Multivariate Regression with Covariance Estimation, MTGS: Multi Traits Genomic Selection.

For any questions: [mbebi@mpimp-golm.mpg.de](mbebi@mpimp-golm.mpg.de) 
