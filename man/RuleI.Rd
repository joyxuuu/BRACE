\name{RuleI}
\alias{RuleI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BRACE RAR Rule I 
}
\description{
Calculate the number of patients allocated to each arm at the interim stage following BRACE RAR Rule I.
}
\usage{
RuleI(nj, n0, pvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nj}{
Total sample size at jth interim stage
}
  \item{n0}{
Sample size for the control arm 
}
  \item{pvec}{
Vector of allocation probability for all treatment arms
}
}
\value{
A vector of number of patients allocated to each arm (for 4 arms scenario only).
}
\examples{
nj <- 200
ctrl_prop <- 1/4
pvec <- c(0.8,0.1,0.1)
n_assign <- RuleI(nj=nj,n0=nj*ctrl_prop,pvec=pvec)		
}
