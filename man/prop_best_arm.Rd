\name{prop_best_arm}
\alias{prop_best_arm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the proportion of patients assigned to the best arm
}
\description{
Calculate the proportion of patients assigned to the best arm in response-adaptive randomization (RAR) design.  
}
\usage{
prop_best_arm(data_all, data_best_arm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_all}{
A list of dataframes from all arms
}
  \item{data_best_arm}{
A dataframe from the best arm
}
}
\value{
The proportion of patients assigned to the best arm, ranging from 0 to 1. 
}
\examples{
\dontrun{
  prop_best_arm(data_all,data_best_arm)
}
}
