prop_best_arm <- function(data_all,data_best_arm){
  return( nrow(data_best_arm)/sum(sapply(data_all,nrow)) )
}

