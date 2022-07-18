\name{mcmc_theta}
\alias{mcmc_theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC function to generate posterior of theta
}
\description{
Run MCMC and generate posterior samples of the estimand (mean transformed OSFD theta)
}
\usage{
mcmc_theta(d, tau, B = 2000, burnin = 500, thin = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
Data input of transformed OSFD 
}
  \item{tau}{
Data input of death indicator
}
  \item{B}{
Number of MCMC iterations
}
  \item{burnin}{
Number of burn-in iterations
}
  \item{thin}{
A value thin = k would result in keeping every kth value and discarding all other values of the MCMC chain.  
}
}
\value{
A vector of theta values of all kept MCMC samples. 
}
\examples{
\dontrun{
theta_ctrl <- mcmc_theta(data_ctrl$d,data_ctrl$tau)
}
}
