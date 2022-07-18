\name{sim_osfd_data}
\alias{sim_osfd_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate organ support-free days (OSFD) data
}
\description{
Simulate organ support-free days (OSFD) data with given sample size and parameter values.
}
\usage{
sim_osfd_data(n, mu, sigma2, omega, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size
}
  \item{mu}{
Mean of transformed OSFD among non-censored patients
}
  \item{sigma2}{
Variance of transformed OSFD among non-censored patients
}
  \item{omega}{
Censoring prportion of patients
}
  \item{lambda}{
Mortality rate 
}
}
\value{
Simulated data consisting of OSFD on original scale (y), death indicator (tau) and transformed OSFD (d).
}
\examples{
n0 <- 50
mean_osfd <- 20
mu_ctrl_true <- -log(30-mean_osfd)
sigma2_true <- 0.8^2
omega_true <- 0.3
lambda_ctrl_true <- 0.2
data_ctrl <- sim_osfd_data(n=n0, mu=mu_ctrl_true, sigma2=sigma2_true, 
omega=omega_true, lambda=lambda_ctrl_true)
}
