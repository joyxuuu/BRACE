mcmc_theta <- function(d,tau,B=2000,burnin=500,thin=10){
  mcmcout <- run_mcmc(d,tau,B=B,burnin=burnin,thin=thin)
  mu <- mcmcout[,"mu"]
  sigma2 <- mcmcout[,"sigma2"]
  sd <- sqrt(sigma2)
  omega <- mcmcout[,"omega"]
  lambda <- mcmcout[,"lambda"]
  rtn <- (dnorm((-log(2)-mu)/sd)-dnorm((-log(30)-mu)/sd))/
    (pnorm((-log(2)-mu)/sd)-pnorm((-log(30)-mu)/sd))
  theta <- (-log(2)*omega +  (1-omega)*(mu-sd*rtn))
  return(theta)
}