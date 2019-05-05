Perform Bayesian Gaussian graphical model with clustering structure for ordered variables in multilevel data sets. The conditional independence between variables are measured on the latent scale via the extended rank likelihood method. Shrinkage effects are applied on the precision matrix to handle multicollinearity. Clustering effects are modelled through the random effects. Missing data are allowed.
The main function 'GCMlasso' fits Bayesian Gaussian copula model with graphical lasso prior for variables with ordering (continuous, ordinal and binary) in multilevel data sets. Adaptive graphical lasso prior is put on the precion matrix of the latent variables conditional on the random effects, where the latent variables are implied by the extended rank likelihood method.
Functions 'compare_group', 'plot_graph', 'reg_coef' and 'predict' perform post-analyses after fitting GCMlasso. 'compare_group' computes the differences of the random effects of variables from two groups. 'plot_graph' plots the conditional dependence between variables on the latent scale where the assocations are computed from the precision matrices from the MCMC. 'reg_coef' computes the regression coefficients of the predictors on a response on the latent variables scale. 'predict' makes predictions for positive outcomes for binary response.
The use of the package is demonstrated by the Framingham heart disease data set, which includes 16 ordinal variables and some missing data. Predictions are made for the binary variable – TenYearCHD(10 year risk of coronary heart disease CHD).


### Installation: 
library(devtools)\
install_github("jialiwang1211/GCMlasso")

### Examples:
### run main function
library(GCMlasso)\
GCMlasso_obj<-GCMlasso(data=Framingham,var_ord=1:15,var_group=16,
  nsamp=1000,odens=1,nwarm=500,seed=1,s=1e-2,t=1e-2,verb=TRUE)
  
### compare variables in cluster 1(education=1&2) and cluster 2(education=3&4)
compare_group(GCMlasso_obj,grp1=1:2,grp2=3:4,var=1:15,credible_level=0.95)
 
### plot conditional dependence between variables
plot_graph(GCMlasso_obj,var=1:15,edge_perc=0.65)

### compute regression coefficients on 'TenYearCHD'
reg_coef(GCMlasso_obj,var_pred=1:14,var_response=15)

### predict probabilities of 'TenYearCHD'
predict_val<-predict(GCMlasso_obj,var_response=15,var_group=16)

