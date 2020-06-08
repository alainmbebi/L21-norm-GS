# L21-norm-GS
R functions to perform ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression and genomic prediction.

L21_featselect.R is a script to run the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression as in [(Nie et all.)]()

CV_L21_featselect.R is a script to run K-folds cross-validation to select the tuning parameter to be used in L21_featselect.R

RidgeEst.R is a script to compute ridge estimate

CV_Ridge.R runs K-folds cross-validation to select the tuning parameter to be used in RidgeEst

L21_joint_estim.R performs the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression and jointly estimate the regression coefficient and precision matrices

CV_L21_joint_estim.R is a script to run K-folds cross-validation to select the tuning parameters to be used in L21_joint_estim.R

MOR.R is a script to run multiple output regression [(He et all. 2016)](doi: 10.1093/bioinformatics/btw249)

CV_MOR is a script to run K-folds cross-validation to select the tuning parameters to be used in MOR.R

By default the all cross-validation scripts uses 5-folds

gwas.r is a script to run GWAS with a mixed model

plots_gwas.r is a script for plotting the results

emma.r use the emma function (Kang et all. 2008) for the MRLE
