# L21-norm-GS

R functions to implement genotype-phenotype association and genomic prediction using ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression, as described in [(A. Mbebi et all. 2020)]().

1- The folder code contains the following R scripts:

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


2- The folder Simulation contains data sets used for simulation. Files are to be readed as in the examples below:

* listX_n50p100s30rho.1 (list of the predictor matrices X for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter \rho=.1)

* listY_n50p100s30rho.1  (list of the response matrices Y for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter \rho=.1)

* trueB_n50p100s30rho.9 (list of the true regression coefficient matrices B for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter \rho=.1)

3- The folder data contains Brassica napus and wheat data sets from [MTGS](https://CRAN.R-project.org/package=MTGS) and [BLRG](https://CRAN.R-project.org/package=BGLR ) packages respectively.

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
