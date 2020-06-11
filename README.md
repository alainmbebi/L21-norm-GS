# L21-norm-GS

R functions to implement genotype-phenotype association and genomic prediction using ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression, as described in [(A. Mbebi et al. 2020)]().

1- The folder code contains the following R scripts:

* L21_featselect.R which uses the model proposed in [(Nie et al. 2010)](http://papers.nips.cc/paper/3988-efficient-and-robust-feature-selection-via-joint-l21-norms-minimization) to implement GS under the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression 

* CV_L21_featselect.R which selects the tuning parameter for L21_featselect.R using K-folds cross-validation

* Ridge_estim.R compute the Ridge estimate

* CV_Ridge.R selects the tuning parameter for Ridge_estim.R using K-folds cross-validation

* L21_joint_estim.R performs GS using the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression that jointly estimates the regression coefficients and precision matrix

* CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

* MOR.R runs multiple output regression [(He et al. 2016)](https://academic.oup.com/bioinformatics/article/32/12/i37/2288681)

* CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

2- The folder Simulation contains data sets used for simulation. Files are to be read as in the examples below:

* listX_n50p100s30rho.1 (list of the predictor matrices X for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter ![equation](https://latex.codecogs.com/gif.latex?%5Crho%3D1))

* listY_n50p100s30rho.1  (list of the response matrices Y for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter ![equation](https://latex.codecogs.com/gif.latex?%5Crho%3D1))

* trueB_n50p100s30rho.9 (list of the true regression coefficient matrices B for all 20 replicates, with n=50, p=100, s=30 and the AR(1) parameter ![equation](https://latex.codecogs.com/gif.latex?%5Crho%3D1))

3- The folder data contains Brassica napus and wheat data sets from [MTGS](https://CRAN.R-project.org/package=MTGS) and [BLRG](https://CRAN.R-project.org/package=BGLR ) packages respectively.

4- The file Example_run.R is an implementation example of all models/functions previously discussed, using the Brassica napus data set. 

5- Notes
* The MTGS.mlasso function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using multivariate LASSO 

* The MTGS.kmlasso function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using Kernelized multivariate LASSO

* The MTGS.mrce function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using MRCE

* By default, all cross-validation scripts use 5-folds

* Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 3.6.1), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.

* The following abbreviations are used, CV: Cross Validation, GBLUP: Genomic Best Linear Unbiased Prediction, GS: Genomic Selection, LASSO: Least Absolute Shrinkage and Selection Operator, MOR: Multiple Output Regression
MRCE: Multivariate Regression with Covariance Estimation, MTGS: Multi Traits Genomic Selection.

6- Licence: GPL-3

For any questions: [mbebi@mpimp-golm.mpg.de](mbebi@mpimp-golm.mpg.de) 
