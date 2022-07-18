\name{success_prob}
\alias{success_prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Success probability of the trial
}
\description{
Calculate the success probability of the trial, defined as the probability of difference in estimand (mean transformed OSFD theta) between the best treatment arm and control arm passing a clinically meaningful difference. 
}
\usage{
success_prob(theta_best, theta_ctrl, clin_diff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta_best}{
The posterior samples of theta in the best treatment arm
}
  \item{theta_ctrl}{
The posterior samples of theta in the control arm
}
  \item{clin_diff}{
Clinically meaningful difference
}
}
\value{
Success probaility of the trial, ranging from 0 to 1. 
}
\examples{
\dontrun{
  clin_diff <-  0.0647669
  success_prob(theta_best,theta_ctrl,clin_diff)
}	
}
