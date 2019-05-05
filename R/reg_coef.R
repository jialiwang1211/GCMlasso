#' Compute regression coefficients
#'
#' Compute the regression coefficients of the predictors on a response on the
#' latent variables scale.
#'
#' @param GCMlasso_obj \code{GCMlasso} object.
#'
#' @param var_pred indices of variables as predictors.
#'
#' @param var_response index of variables as response.
#'
#' @param credible_level level of credible interval.
#'
#' @return Posterior means and credible intervals of the coefficients.
#'
#' @examples reg_coef(GCMlasso_obj,var_pred=1:14,var_response=15)
#'
#' @export
reg_coef<-function(GCMlasso_obj,var_pred,var_response,credible_level=0.95){
  data<-GCMlasso_obj$data_ordered
  C<-GCMlasso_obj$Gamma.st
  numsamp=dim(C)[3]
  coef<-matrix(0,length(var_pred),numsamp)

  for (i in 1:numsamp){
    coef[,i]<-(C[var_response,var_pred,i])%*%solve(C[var_pred,var_pred,i])
  }

  mean<-apply(coef,1,mean)
  CI_L<-CI_R<-rep(0,length(var_pred))
  for (i in 1:length(var_pred)){
    CI_L[i]<-HPDinterval(as.mcmc(coef[i,]), prob=credible_level)[1]
    CI_R[i]<-HPDinterval(as.mcmc(coef[i,]), prob=credible_level)[2]
  }

  table<-data.frame(variable=names(data)[var_pred],mean=round(mean,3),
                    CI=paste("(",round(CI_L,3),",",round(CI_R,3),")",sep=""))
  return(table)
}
