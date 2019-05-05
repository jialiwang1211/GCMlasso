#' Gaussian graphical model for ordered variables
#'
#' Perform Bayesian Gaussian graphical model for ordered variables
#' with clustering structure. The conditional independence between variables
#' are measured on the latent scale via the extended rank likelihood method.
#' Shrinkage effects are applied on the precision matrix to handle
#' multicollinearity. Clustering effects are modelled through the random effects.
#' Missing data are allowed.
#'
#' \code{GCMlasso} function fits the Bayesian Gaussian copula model with graphical lasso prior
#' for variables with ordering (continuous, ordinal and binary) in multilevel data sets.
#' Adaptive graphical lasso prior is put on the precion matrix of the latent variables
#' conditional on the random effects, where the latent variables are implied by the
#' extended rank likelihood method.
#'
#' \code{var_group} is the index of the varible that defines the clusters, and
#' should be placed at the last column in the \code{data}. The coding for the clustering variable
#' is from 1 to the total number of clusters. The binary variables
#' in \code{var_ord} should be coded as 0 (control) or 1(case) in the case-control studies.
#'
#' Missing data are allowed for ordered variables and should be denoted as NA.
#'
#' @param data an N by p data frame,
#'
#' @param var_ord indices of ordinal variables.
#'
#' @param var_group index of the variable that defines clusters.
#'
#' @param nsamp number of MCMC iterations.
#'
#' @param odens number of iterations between saved samples.
#'
#' @param nwarm number of MCMC iterations as burn-in.
#'
#' @param seed a random integer.
#'
#' @param s hyperparameter of \eqn{lambda}, degree of freedom in adaptive graphical lasso prior.
#'
#' @param t hyperparameter of \eqn{lambda}, shrinkage in adaptive graphical lasso prior.
#'
#' @param verb print progress of MCMC, logical TRUE or FALSE.
#'
#' @return An object with S3 class "\code{GCMlasso}" is returned.
#'
#' \item{data_ordered}{the same as the input data but ordered by the
#' clustering variable.}
#'
#' \item{Gamma.st}{saved variance covariance matrices for latent variables.}
#'
#' \item{Omega.st}{saved precision matrices for latent variables.}
#'
#' \item{psi.st}{saved variance covariance matrices for random effects.}
#'
#' \item{z.st}{saved latent variables.}
#'
#' \item{b.st}{saved random effects.}
#'
#' @examples GCMlasso_obj<-GCMlasso(data=Framingham,var_ord=1:15,var_group=16,
#'   nsamp=1000,odens=1,nwarm=500,seed=1,s=1e-2,t=1e-2,verb=TRUE)
#' @author Jiali Wang (\email{jiali.wang@@data61.csiro.au})
#'
#' @references
#' \insertRef{hoff2007extending}{GCMlasso}
#'
#' \insertRef{wang2012bayesian}{GCMlasso}
#'
#' \insertRef{wang2017copula}{GCMlasso}
#'
#' @export
#' @import MASS
#' @import MCMCpack
#' @importFrom Rdpack reprompt
GCMlasso<-function(data,var_ord,var_group,nsamp=1000,odens=1,nwarm=500,
                      seed=1,s=1e-2,t=1e-2,verb=TRUE){
  set.seed(seed)

  data <- data[order(data[,var_group]),]
  Y<-data[,var_ord]
  ni<-as.numeric(table(data[,var_group]))
  m<- length(ni)
  N<-nrow(Y)
  group<-rep(1:m,ni)
  p<-dim(Y)[2]

  nu_psi<-p+2
  lambda_psi_inv<-solve(100*diag(p))

  lambda_ii<-1

  R<-NULL
  for(k in 1:p) {
    R<-cbind(R, match(Y[,k],sort(unique(Y[,k]))))
  }
  Rlevels<-apply(R,2,max,na.rm=TRUE)
  Ranks<-apply(Y[,1:p],2,rank,ties.method="random",na.last="keep")

  N.full<-apply(!is.na(Ranks),2,sum)
  U<-t(t(Ranks)/(N.full+1))
  Z<-qnorm(U)
  Zfill<-matrix(rnorm(N*p),N,p)
  Z[is.na(Y[,1:p])]<-Zfill[is.na(Y[,1:p])]

  C.st<-Omega.st<-psi.st<-array(dim=c(p,p,floor((nsamp-nwarm)/odens)))
  b.st<-array(dim=c(m,p,floor((nsamp-nwarm)/odens)))
  z.st<- array(dim=c(N,p,floor((nsamp-nwarm)/odens)))

  C<-Omega<-psi<-diag(p)
  tau <- matrix(0,p,p)

  indmx<-matrix(1:p^2,p,p)
  upperind<-indmx[upper.tri(indmx)]
  indmx_t <- t(indmx)
  lowerind <- indmx_t[upper.tri(indmx_t)]

  ind_noi_all<-matrix(0,p-1,p)
  for (i in 1:p){
    ind_noi_all[,i]<-(1:p)[-i]
  }


  for (ns in 1:nsamp){

    b<-NULL
    for (i in 1:m){
      vb<-solve(ni[i]*solve(C)+solve(psi))
      mb<-vb%*%solve(C)%*%apply(Z[group==i,],2,sum)
      b<-rbind(b,mvrnorm(1,mb,vb))
    }

    b.expand<-NULL
    for (i in 1:m){
      b.expand<- rbind(b.expand,kronecker(t(b[i,]),rep(1,ni[i])))
    }

    psi<-riwish(nu_psi+m,lambda_psi_inv+t(b)%*%b)

    S<-t(Z-b.expand)%*%(Z-b.expand)

    Omegaadjust <- pmax(abs(Omega[upperind]),1e-12)
    lambda<-rgamma(length(lowerind),shape=1+s,scale=1/(Omegaadjust+t))

    lambda_prime<-lambda^2
    mu_prime<-pmin(lambda/Omegaadjust,1e12)
    tau_temp<-pmax(1/rinvgauss(length(lowerind),mu_prime,lambda_prime),1e-12)
    tau[upperind]<-tau_temp
    tau[lowerind]<-tau_temp

    for (k in 1:p){
      ind_noi <- ind_noi_all[,k]
      tau_temp <- tau[ind_noi,k]
      invOmega11 <-solve(Omega[ind_noi,ind_noi])

      Omegai <- (S[k,k]+lambda_ii)*invOmega11 +diag(1/tau_temp)
      beta<-mvrnorm(1,-solve(Omegai)%*%S[ind_noi,k], solve(Omegai))

      Omega[ind_noi,k]<-beta
      Omega[k,ind_noi]<-beta
      gam<-rgamma(1,shape=N/2+1,scale=2/(S[k,k]+lambda_ii))

      Omega[k,k]<-gam+t(beta)%*%invOmega11%*%beta
    }

    C<-solve(Omega)
    diag<-sqrt(diag(C+psi))
    C<-t((C)/diag)/diag
    psi<-t((psi)/diag)/diag

    for (k in 1:p) {
      Skc<- C[k,-k]%*%solve(C[-k,-k])
      sdk<- sqrt(C[k,k] -C[k,-k]%*%solve(C[-k,-k])%*%C[-k,k])
      muk<- (Z[,-k]-b.expand[,-k])%*%t(Skc)+b.expand[,k]

      for(r in 1:Rlevels[k]){
        ir<- (1:N)[R[,k]==r & !is.na(R[,k])]
        lb<-suppressWarnings(max( Z[ R[,k]==r-1,k],na.rm=TRUE))
        ub<-suppressWarnings(min( Z[ R[,k]==r+1,k],na.rm=TRUE))
        Z[ir,k]<-qnorm(runif(length(ir),pnorm(lb,muk[ir],sdk),pnorm(ub,muk[ir],sdk)),muk[ir],sdk)
      }

      ir<-(1:N)[is.na(R[,k])]
      Z[ir,k]<-rnorm(length(ir),muk[ir],sdk)
    }

    if (((ns-nwarm)%%odens==0)&(ns>=nwarm)){
      C.st[,,(ns-nwarm)/odens]<-C
      Omega.st[,,(ns-nwarm)/odens]<-Omega
      b.st[,,(ns-nwarm)/odens]<-b
      psi.st[,,(ns-nwarm)/odens]<-psi
      z.st[,,(ns-nwarm)/odens]<-Z
    }

    if (verb == TRUE & (ns%%(odens * 100)) == 0) {
      cat(round(100 * ns/nsamp), "% done \n")
    }
  }
  output<-list(data_ordered=data,Gamma.st=C.st, Omega.st=Omega.st, psi.st=psi.st,z.st=z.st,b.st=b.st)
}

