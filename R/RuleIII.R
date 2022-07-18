RuleIII <- function(nj,pvec){
    
    pnorm <- pvec/sum(pvec)
    pbest <- max(pnorm)
    scale_factor <- 1/pbest
    
    p0 <- pbest 
    pvec2 <- c(p0,pnorm)
    pvec2[-1][-which.max(pnorm)] <- scale_factor*pnorm[-which.max(pnorm)]
    pnorm2 <- pvec2/sum(pvec2)
    
    nvec = round(nj*pnorm2)
    
    return(nvec)

}    