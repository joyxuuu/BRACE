RuleII <- function(nj,pvec){
  
  cap <- round(nj*0.8)
  
    nvec <- rep(NA,3)
    pnorm <- pvec/sum(pvec)
    
    pbest <- max(pnorm)
    nvec[which.max(pnorm)] <- min(round(nj*pbest),cap)
    
    n0 = round((nj - nvec[which.max(pnorm)])*(1/(length(pvec))))
    
    pvec_sub <- pnorm[-which.max(pnorm)]
    if(sum(pvec_sub)==0){
      nvec[-which.max(pnorm)] <- 0
    } else{
      nvec[-which.max(pnorm)] <- round((nj - nvec[which.max(pnorm)]-n0)*(pvec_sub/sum(pvec_sub)))
    }
    nvec2 <- c(n0,nvec)  
    return(nvec2)

}    
  