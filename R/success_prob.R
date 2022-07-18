success_prob <- function(theta_best,theta_ctrl,clin_diff){
  numerator <- outer(theta_best,theta_ctrl,FUN=function(x,y) x-y>clin_diff)
  return(sum(numerator)/length(numerator))
}