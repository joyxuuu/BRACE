RuleI <- function(nj,n0,pvec){
    nvec = round((nj-n0)*pvec/sum(pvec))
    return(c(n0,nvec))
}