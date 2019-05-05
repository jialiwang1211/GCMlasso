#' Sample from inverse Gaussian distribution
#'
#' Draw samples from inverse Gaussian distribution
#'
#' @param n sample size.
#' @param mean vector of positive means.
#' @param shape vector of (positive) shape parameters.
#' @param dispersion vector of (positive) dispersion parameters.
#'
#' @return A vector of random samples from inverse Gaussian distribution
#' @examples rinvgauss(100,mean=1,shape=2)
#' @export
rinvgauss<-function (n, mean = 1, shape = NULL, dispersion = 1)
{
  if (!is.null(shape))
    dispersion <- 1/shape
  if (length(n) > 1L)
    n <- length(n)
  else n <- as.integer(n)
  if (n < 0L)
    stop("n can't be negative")
  if (n == 0L || length(mean) == 0L || length(dispersion) ==
      0L)
    return(numeric(0L))
  mu <- rep_len(mean, n)
  phi <- rep_len(dispersion, n)
  r <- rep_len(0, n)
  mu.ok <- (mu > 0 & is.finite(mu))
  phi.ok <- (phi > 0 & is.finite(phi))
  i <- (mu.ok & phi.ok)
  if (!all(i)) {
    j <- !i
    invchisq <- (mu[j] == Inf & phi.ok[j])
    invchisq[is.na(invchisq)] <- FALSE
    if (any(invchisq)) {
      m <- sum(invchisq)
      r[j][invchisq] <- rnorm(m)^(-2)/phi[j][invchisq]
      j[j][invchisq] <- FALSE
    }
    infdisp <- (phi[j] == Inf)
    infdisp[is.na(infdisp)] <- FALSE
    if (any(infdisp)) {
      r[j][infdisp] <- 0
      j[j][infdisp] <- FALSE
    }
    r[j] <- NA
    n <- sum(i)
    if (n == 0L)
      return(r)
  }
  Y <- rnorm(n)^2
  Yphi <- Y * phi[i] * mu[i]
  bigphi <- (Yphi > 5e+05)
  if (any(bigphi)) {
    X1 <- Y
    X1[bigphi] <- 1/Yphi[bigphi]
    X1[!bigphi] <- 1 + Yphi[!bigphi]/2 * (1 - sqrt(1 + 4/Yphi[!bigphi]))
  }
  else {
    X1 <- 1 + Yphi/2 * (1 - sqrt(1 + 4/Yphi))
  }
  firstroot <- (runif(n) < 1/(1 + X1))
  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  r[i] <- mu[i] * r[i]

  return(r)
}
