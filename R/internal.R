#############################
###### Internal functions ###
#############################

### run MCMC 

### constants:

alpha0 <- 1e-4
beta0 <- 1e-4
mu0 <- 0
sigma2_0 <- 10^4
lb <- 30 ## -log(30) is LB
ub <- 2 ## -log(2) is UB

### updating functions 

update_ldnt <- function(ld,mu,sigma2,eps=1e-5){
  ldstd <- (ld - mu)/sqrt(sigma2)
  x <- pnorm(ldstd)
  a <- pnorm((-log(lb)-mu)/sqrt(sigma2))
  b <- pnorm((-log(ub)-mu)/sqrt(sigma2))
  
  ldnt <- mu + sqrt(sigma2)*qnorm((x-a+eps)/(b-a+eps)) 
  return(ldnt)
}

update_mu <- function(ld,sigma2){
  n <- length(ld)
  mean = (1/(1/sigma2_0 + n/sigma2))*(mu0/sigma2_0 + sum(ld)/sigma2)
  sd = sqrt(1/(1/sigma2_0 + n/sigma2)) 
  return(rnorm(1,mean=mean,sd=sd))
}

update_sigma2 <- function(ld,mu){
  n <- length(ld)
  shape=alpha0+n/2
  scale=beta0+sum((ld-mu)^2)/2
  return(rinvgamma(1, shape=shape, scale = scale))
}

update_omega <- function(z){
  shape1=1+sum(z)
  shape2=1+length(z)-sum(z)
  return(rbeta(1,shape1,shape2))
}

update_lambda <- function(tau){
  shape1=1+sum(tau)
  shape2=1+length(tau)-sum(tau)
  return(rbeta(1,shape1,shape2))
}


run_mcmc <- function(d,tau,B,burnin,thin){
  
  ### constants
  alpha0 <- 1e-4
  beta0 <- 1e-4
  mu0 <- 0
  sigma2_0 <- 10^4
  lb <- 30 ## -log(30) is LB
  ub <- 2 ## -log(2) is UB
  

  n <- length(tau)
  d_live <- d[tau==0]
  n_live <- sum(tau==0)
  
  ind_d_tn <- which(d_live <= -log(3)) 
  ind_d_peak <-   which(d_live > -log(3))
  
  d_tn <- d_live[ind_d_tn]

  mu.init <- mean(d_tn)
  sigma2.init <- var(d_tn)
  omega.init <- 0.5
  lambda.init <- (n-n_live)/n
  
  store.mat <- matrix(NA,nrow=B,ncol=4)
  colnames(store.mat) <- c("mu","sigma2","omega","lambda")
  
  mu <- store.mat[1,"mu"] <- mu.init
  sigma2 <- store.mat[1,"sigma2"] <- sigma2.init
  omega <- store.mat[1,"omega"] <- omega.init
  lambda <- store.mat[1,"lambda"] <- lambda.init
  
  z <- rbinom(n_live,1,omega)
  ind_tn <- which(z==0)
  
  for(b in 2:B){
    mu <- store.mat[b,"mu"] <- update_mu(ld=d_live[ind_tn],sigma2=sigma2)
    sigma2 <- store.mat[b,"sigma2"] <- update_sigma2(ld=d_live[ind_tn],mu=mu)
    
    z <- rep(0,n_live)
    zz <- try(rbinom(length(ind_d_peak),1,
                     prob=omega/(omega+(1-omega)*dtruncnorm(x=-log(2),a=-log(lb),b=-log(ub), mean=mu,sd=sqrt(sigma2)))),silent=T)
    
    if(class(zz)== "try-error"){
      z[ind_d_peak] <- rbinom(length(ind_d_peak),1,omega)
    } else{
      z[ind_d_peak] <- zz
    }
    
    ind_tn <- which(z==0)
    
    omega <- store.mat[b,"omega"] <- update_omega(z=z)
    lambda <- store.mat[b,"lambda"] <- update_lambda(tau=tau)
    
  }
  
  store.mat2 <- store.mat[seq(burnin+1,B,thin),]
  
  return(store.mat2)
  
}



