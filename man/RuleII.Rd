\name{RuleII}
\alias{RuleII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BRACE RAR Rule II 
}
\description{
Calculate the number of patients allocated to each arm at the interim stage following BRACE RAR Rule II.
}
\usage{
RuleII(nj, pvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nj}{
Total sample size at jth interim stage
}
  \item{pvec}{
Vector of allocation probability for all treatment arms
}
}
\value{
A vector of number of patients allocated to each arm.
}
\examples{
nj <- 200
pvec <- c(0.8,0.1,0.1)
n_assign <- RuleII(nj=nj,pvec=pvec)	
}
