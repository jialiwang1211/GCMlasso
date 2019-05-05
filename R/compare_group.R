#' Compare variables in two groups
#'
#' Compute the differences of the random effects for variables from two groups.
#'
#' @param GCMlasso_obj \code{GCMlasso} object.
#'
#' @param grp1 levels in the cluster variable that belong to group1.
#'
#' @param grp2 levels in the cluster variable that belong to group2.
#'
#' @param var variables to be compared.
#'
#' @param credible_level level of the credible interval.
#'
#' @return Posterior means and credible intervals of the differences
#' of variables in two groups.
#' @examples compare_group(GCMlasso_obj=GCMlasso_obj,grp1=1:2,grp2=3:4,
#'  var=1:15,credible_level=0.95)
#' @export
#' @import coda
compare_group<-function(GCMlasso_obj=GCMlasso_obj,grp1,grp2,var,credible_level=0.95){
  data<-GCMlasso_obj$data_ordered
  b<-GCMlasso_obj$b.st[,var,]
  numsamp=dim(b)[3]
  b.grp1<-b.grp2<-matrix(0,length(var),numsamp)
  for (i in 1:numsamp){
    for (j in 1:length(var)){
      b.grp1[j,i]<-mean(b[grp1,j,i])
      b.grp2[j,i]<-mean(b[grp2,j,i])
    }}
  CI_L<-CI_R<-mean<-rep(0,length(var))
  mean<-apply(b.grp1-b.grp2,1,mean)
  for (j in 1:length(var)){
    CI_L[j]<-HPDinterval(as.mcmc(b.grp1[j,]-b.grp2[j,]), prob=credible_level)[1]
    CI_R[j]<-HPDinterval(as.mcmc(b.grp1[j,]-b.grp2[j,]), prob=credible_level)[2]
  }
  table<-data.frame(variable=names(data)[var],mean=round(mean,3),
                    CI=paste("(",round(CI_L,3),",",round(CI_R,3),")",sep=""))
  return(table)
}
