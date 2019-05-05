#' Predict positive outcomes
#'
#' Making predictions for positive outcomes for binary response.
#'
#' Include new observations in the data set by treating the predicting variables
#' as NA, and run the MCMC using the function \code{GCMlasso}. Code positive
#' cases as 1 and controls as 0. The maximum latent variable corresponding to the
#' positive cases is used as threshold to predict new observations.
#'
#' @param GCMlasso_obj \code{GCMlasso} object.
#'
#' @param var_response index of binary response variable.
#'
#' @param var_group index of the variable that defines clusters.
#'
#' @param rep_sample number of samples in the posterior predictive sampling.
#'
#' @param seed a random integer.
#'
#' @return Mean values of posterior predictive probabilities of positive cases
#' for observations with known outcomes and new observations.
#' @examples predict_val<-predict(GCMlasso_obj,var_response=15,var_group=16)
#' @export
predict<-function(GCMlasso_obj,var_response,var_group,rep_sample=100,seed=1){
  set.seed(seed)
  data<-GCMlasso_obj$data_ordered
  z<-GCMlasso_obj$z.st
  b<-GCMlasso_obj$b.st
  C<-GCMlasso_obj$Gamma.st

  numsamp<-dim(z)[3]
  ni<-as.numeric(table(data[,var_group]))
  m<- length(ni)

  rep_pred<-matrix(0,nrow(data),rep_sample)
  idx<-sample(1:numsamp,rep_sample)

  for (rep in 1:rep_sample){
    id<-idx[rep]
    b.expand<-NULL
    for (i in 1:m){
      b.expand<- rbind(b.expand,kronecker(t(b[i,,id]),rep(1,ni[i])))
    }

    sdk<- sqrt(C[var_response,var_response,id]-C[var_response,-var_response,id]
               %*%solve(C[-var_response,-var_response,id])%*%C[-var_response,var_response,id])

    z_sample<-y_pred<-rep(0,nrow(data))

    for (i in 1:nrow(data)){
      z_sample[i]<-mvrnorm(1,(z[i,-var_response,id]-b.expand[i,-var_response])
                           %*%t(C[var_response,-var_response,id]%*%solve(C[-var_response,-var_response,id]))+
                             b.expand[i,var_response],sdk)
    }

    z_save<-z[,var_response,id]
    y_pred[z_sample>min(z_save[data[,var_response]==1],na.rm = TRUE)]<-1

    rep_pred[,rep]<-y_pred
  }

  predscore<-apply(rep_pred,1,mean)
  return(predscore=predscore)
}
