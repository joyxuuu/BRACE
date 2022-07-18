sim_osfd_data <- function(n,mu,sigma2,omega,lambda){
  
  #generate tau 
  tau <- rbinom(n,1,lambda)
  n_live <- sum(tau==0)
  ind_live <- which(tau==0)
  
  #generate d and y 

  bin_mix <- rbinom(n_live,1,prob=omega) ## 0=from TN component, 1=from censored component
  n_nonone <- sum(bin_mix==0)
  n_one <- n_live - n_nonone

  if(n_one==0) {
    d <- y <- rep(0,n)
    d[ind_live] <- rtruncnorm(n_nonone,a=-log(30),b=-log(2),
                                mean = mu, 
                                sd = sqrt(sigma2))
    y[ind_live] <- round(30-exp(-d[ind_live]))
  } else if (n_nonone==0){
    d <- y <- rep(0,n)
    d[ind_live] <- -log(2)
    y[ind_live] <- 28
  } else {
    d <- y <- rep(0,n)
    d[ind_live][which(bin_mix==0)] <- rtruncnorm(n_nonone,a=-log(30),b=-log(2),
                                                    mean = mu, sd = sqrt(sigma2))
    d[ind_live][which(bin_mix==1)] <- -log(2) 
    y[ind_live][which(bin_mix==0)] <- round(30-exp(-d[ind_live][which(bin_mix==0)])) 
    y[ind_live][which(bin_mix==1)] <- 28
    
  }

  data <- data.frame(y=y,tau=tau,d=d)
  
  return(data)
}

