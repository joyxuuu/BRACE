## test BRACE ###

library(BRACE)

## examples ##

# True parameter setting:

n_total <- 2000
K <- 4 # arm k=1 is control arm, arm k=2-4 are treatment
J <- 10 
nj <- n_total/J # n=200 in each stage 
ctrl_prop <- 1/4


lambda_ctrl_true <- 0.2
lambda_trt1_true <- 0.15
lambda_trt2_true <- 0.2
lambda_trt3_true <- 0.2

omega_true <- 0.3

mean_osfd <- 20
mu_ctrl_true <- -log(30-mean_osfd)
es <- 0.8 
sigma2_true <- 0.8^2
mu_trt1_true <- mu_ctrl_true + (es*sqrt(sigma2_true))
mu_trt2_true <- mu_ctrl_true 
mu_trt3_true <- mu_ctrl_true 

R <- 10 # number of replications

final_theta_ctrl_R1 <- final_theta_ctrl_R2 <- final_theta_ctrl_R3 <- 
final_theta_trt1_R1 <- final_theta_trt1_R2 <- final_theta_trt1_R3 <- 
final_theta_trt2_R1 <- final_theta_trt2_R2 <- final_theta_trt2_R3 <- 
final_theta_trt3_R1 <- final_theta_trt3_R2 <- final_theta_trt3_R3 <- vector("list",length=R)

final_data_ctrl_R1 <- final_data_ctrl_R2 <- final_data_ctrl_R3 <- 
final_data_trt1_R1 <- final_data_trt1_R2 <- final_data_trt1_R3 <-  
final_data_trt2_R1 <- final_data_trt2_R2 <- final_data_trt2_R3 <- 
final_data_trt3_R1 <- final_data_trt3_R2 <- final_data_trt3_R3 <- vector("list",length=R)


##### calculate allocation prob

prob_best_3trtarm <- function(theta1,theta2,theta3){
  # by default assesses: P(theta1<theta2 & theta1<theta3)
  numerator <- as.numeric(outer(theta1,theta2,">"))*as.numeric(outer(theta1,theta3,">"))
  return(sum(numerator)/length(numerator))
}

set.seed(12345)

for(r in 1:R){
  
  print(paste("Simulation",r,sep=" "))
  start <- proc.time()
  
  ## Rule I
  
  print("Rule I")
  
  #initial stage 
  
  n0 <- n1 <- n2 <- n3 <- nj*(1/4) ## initial ss
  
  data_init_ctrl <- data_ctrl <- sim_osfd_data(n=n0,mu=mu_ctrl_true,sigma2_true,
                                             omega_true,lambda_ctrl_true)
  data_init_trt1 <- data_trt1 <- sim_osfd_data(n=n1,mu=mu_trt1_true,sigma2_true,
                                             omega_true,lambda_trt1_true)
  data_init_trt2 <- data_trt2 <- sim_osfd_data(n=n2,mu=mu_trt2_true,sigma2_true,
                                             omega_true,lambda_trt2_true)
  data_init_trt3 <- data_trt3 <- sim_osfd_data(n=n3,mu=mu_trt3_true,sigma2_true,
                                             omega_true,lambda_trt3_true)
  
  theta_ctrl <- mcmc_theta(data_ctrl$d,data_ctrl$tau)
  theta_trt1 <- mcmc_theta(data_trt1$d,data_trt1$tau)
  theta_trt2 <- mcmc_theta(data_trt2$d,data_trt2$tau)
  theta_trt3 <- mcmc_theta(data_trt3$d,data_trt3$tau)
 
  p1_init <- p1 <- prob_best_3trtarm(theta_trt1,theta_trt2,theta_trt3)
  p2_init <- p2 <- prob_best_3trtarm(theta_trt2,theta_trt1,theta_trt3)
  p3_init <- p3 <- prob_best_3trtarm(theta_trt3,theta_trt1,theta_trt2)
  
  ## interim stage
  
  for(j in 2:J){
    
    n_assign <- RuleI(nj=nj,n0=nj*ctrl_prop,pvec=c(p1,p2,p3))
  
    n0 <- n_assign[1]
    n1 <- n_assign[2]
    n2 <- n_assign[3]
    n3 <- n_assign[4]

    if(n0!=0){
      dataj_ctrl <- sim_osfd_data(n=n0,mu=mu_ctrl_true,sigma2_true,
                                  omega_true,lambda_ctrl_true)
    } else{
      dataj_ctrl <- NULL }
    
    if(n1!=0){
      dataj_trt1 <- sim_osfd_data(n=n1,mu=mu_trt1_true,sigma2_true,
                                  omega_true,lambda_trt1_true)
    } else{
      dataj_trt1 <- NULL } 
    
    if(n2!=0){
      dataj_trt2 <- sim_osfd_data(n=n2,mu=mu_trt2_true,sigma2_true,
                                  omega_true,lambda_trt2_true)
    } else{
      dataj_trt2 <- NULL } 
    
    if(n3!=0){
      dataj_trt3 <- sim_osfd_data(n=n3,mu=mu_trt3_true,sigma2_true,
                                  omega_true,lambda_trt3_true)
    } else{
      dataj_trt3 <- NULL }
    
    if(!is.null(dataj_ctrl)){
        data_ctrl <- rbind(data_ctrl,dataj_ctrl)}
    
    if(!is.null(dataj_trt1)){
        data_trt1 <- rbind(data_trt1,dataj_trt1)}
    
    if(!is.null(dataj_trt2)){
        data_trt2 <- rbind(data_trt2,dataj_trt2)}

    if(!is.null(dataj_trt3)){
        data_trt3 <- rbind(data_trt3,dataj_trt3)}
    
    
    theta_ctrl <- mcmc_theta(data_ctrl$d,data_ctrl$tau)
    theta_trt1 <- mcmc_theta(data_trt1$d,data_trt1$tau)
    theta_trt2 <- mcmc_theta(data_trt2$d,data_trt2$tau)
    theta_trt3 <- mcmc_theta(data_trt3$d,data_trt3$tau)
    
    p1 <- prob_best_3trtarm(theta_trt1,theta_trt2,theta_trt3)
    p2 <- prob_best_3trtarm(theta_trt2,theta_trt1,theta_trt3)
    p3 <- prob_best_3trtarm(theta_trt3,theta_trt1,theta_trt2)
    
    n0 <- n1 <- n2 <- n3 <- nj*ctrl_prop 
    
    rm(dataj_ctrl)
    rm(dataj_trt1)
    rm(dataj_trt2)
    rm(dataj_trt3)
    
  }
  
  final_data_ctrl_R1[[r]] <- data_ctrl
  final_data_trt1_R1[[r]] <- data_trt1
  final_data_trt2_R1[[r]] <- data_trt2
  final_data_trt3_R1[[r]] <- data_trt3
  
  final_theta_ctrl_R1[[r]] <- theta_ctrl
  final_theta_trt1_R1[[r]] <- theta_trt1
  final_theta_trt2_R1[[r]] <- theta_trt2
  final_theta_trt3_R1[[r]] <- theta_trt3
  
  
  ## Rule II
  
  print("Rule II")
  
  
  data_ctrl <- data_init_ctrl
  data_trt1 <- data_init_ctrl
  data_trt2 <- data_init_ctrl
  data_trt3 <- data_init_ctrl
  
  p1 <- p1_init
  p2 <- p2_init
  p3 <- p3_init
  
  for(j in 2:J){
    
    n_assign <- RuleII(nj=nj,pvec=c(p1,p2,p3))
    
    n0 <- n_assign[1]
    n1 <- n_assign[2]
    n2 <- n_assign[3]
    n3 <- n_assign[4]
    
    if(n0!=0){
      dataj_ctrl <- sim_osfd_data(n=n0,mu=mu_ctrl_true,sigma2_true,
                                  omega_true,lambda_ctrl_true)
    } else{
      dataj_ctrl <- NULL }
    
    if(n1!=0){
      dataj_trt1 <- sim_osfd_data(n=n1,mu=mu_trt1_true,sigma2_true,
                                  omega_true,lambda_trt1_true)
    } else{
      dataj_trt1 <- NULL } 
    
    if(n2!=0){
      dataj_trt2 <- sim_osfd_data(n=n2,mu=mu_trt2_true,sigma2_true,
                                  omega_true,lambda_trt2_true)
    } else{
      dataj_trt2 <- NULL } 
    
    if(n3!=0){
      dataj_trt3 <- sim_osfd_data(n=n3,mu=mu_trt3_true,sigma2_true,
                                  omega_true,lambda_trt3_true)
    } else{
      dataj_trt3 <- NULL }
    
    if(!is.null(dataj_ctrl)){
      data_ctrl <- rbind(data_ctrl,dataj_ctrl)}
    
    if(!is.null(dataj_trt1)){
      data_trt1 <- rbind(data_trt1,dataj_trt1)}
    
    if(!is.null(dataj_trt2)){
      data_trt2 <- rbind(data_trt2,dataj_trt2)}
    
    if(!is.null(dataj_trt3)){
      data_trt3 <- rbind(data_trt3,dataj_trt3)}
    
    
    theta_ctrl <- mcmc_theta(data_ctrl$d,data_ctrl$tau)
    theta_trt1 <- mcmc_theta(data_trt1$d,data_trt1$tau)
    theta_trt2 <- mcmc_theta(data_trt2$d,data_trt2$tau)
    theta_trt3 <- mcmc_theta(data_trt3$d,data_trt3$tau)
    
    p1 <- prob_best_3trtarm(theta_trt1,theta_trt2,theta_trt3)
    p2 <- prob_best_3trtarm(theta_trt2,theta_trt1,theta_trt3)
    p3 <- prob_best_3trtarm(theta_trt3,theta_trt1,theta_trt2)
    
    n0 <- n1 <- n2 <- n3 <- nj*ctrl_prop 
    
    rm(dataj_ctrl)
    rm(dataj_trt1)
    rm(dataj_trt2)
    rm(dataj_trt3)
    
  }
  
  final_data_ctrl_R2[[r]] <- data_ctrl
  final_data_trt1_R2[[r]] <- data_trt1
  final_data_trt2_R2[[r]] <- data_trt2
  final_data_trt3_R2[[r]] <- data_trt3

  final_theta_ctrl_R2[[r]] <- theta_ctrl
  final_theta_trt1_R2[[r]] <- theta_trt1
  final_theta_trt2_R2[[r]] <- theta_trt2
  final_theta_trt3_R2[[r]] <- theta_trt3
  
    
  ## Rule III
  
  print("Rule III")
  
  
  data_ctrl <- data_init_ctrl
  data_trt1 <- data_init_ctrl
  data_trt2 <- data_init_ctrl
  data_trt3 <- data_init_ctrl
  
  p1 <- p1_init
  p2 <- p2_init
  p3 <- p3_init
  
  for(j in 2:J){
    
    n_assign <- RuleIII(nj=nj,pvec=c(p1,p2,p3))
    
    n0 <- n_assign[1]
    n1 <- n_assign[2]
    n2 <- n_assign[3]
    n3 <- n_assign[4]
    
    if(n0!=0){
      dataj_ctrl <- sim_osfd_data(n=n0,mu=mu_ctrl_true,sigma2_true,
                                  omega_true,lambda_ctrl_true)
    } else{
      dataj_ctrl <- NULL }
    
    if(n1!=0){
      dataj_trt1 <- sim_osfd_data(n=n1,mu=mu_trt1_true,sigma2_true,
                                  omega_true,lambda_trt1_true)
    } else{
      dataj_trt1 <- NULL } 
    
    if(n2!=0){
      dataj_trt2 <- sim_osfd_data(n=n2,mu=mu_trt2_true,sigma2_true,
                                  omega_true,lambda_trt2_true)
    } else{
      dataj_trt2 <- NULL } 
    
    if(n3!=0){
      dataj_trt3 <- sim_osfd_data(n=n3,mu=mu_trt3_true,sigma2_true,
                                  omega_true,lambda_trt3_true)
    } else{
      dataj_trt3 <- NULL }
    
    if(!is.null(dataj_ctrl)){
      data_ctrl <- rbind(data_ctrl,dataj_ctrl)}
    
    if(!is.null(dataj_trt1)){
      data_trt1 <- rbind(data_trt1,dataj_trt1)}
    
    if(!is.null(dataj_trt2)){
      data_trt2 <- rbind(data_trt2,dataj_trt2)}
    
    if(!is.null(dataj_trt3)){
      data_trt3 <- rbind(data_trt3,dataj_trt3)}
    
    
    theta_ctrl <- mcmc_theta(data_ctrl$d,data_ctrl$tau)
    theta_trt1 <- mcmc_theta(data_trt1$d,data_trt1$tau)
    theta_trt2 <- mcmc_theta(data_trt2$d,data_trt2$tau)
    theta_trt3 <- mcmc_theta(data_trt3$d,data_trt3$tau)
    
    p1 <- prob_best_3trtarm(theta_trt1,theta_trt2,theta_trt3)
    p2 <- prob_best_3trtarm(theta_trt2,theta_trt1,theta_trt3)
    p3 <- prob_best_3trtarm(theta_trt3,theta_trt1,theta_trt2)
    
    n0 <- n1 <- n2 <- n3 <- nj*ctrl_prop 
    
    rm(dataj_ctrl)
    rm(dataj_trt1)
    rm(dataj_trt2)
    rm(dataj_trt3)
    
  }
  
  final_data_ctrl_R3[[r]] <- data_ctrl
  final_data_trt1_R3[[r]] <- data_trt1
  final_data_trt2_R3[[r]] <- data_trt2
  final_data_trt3_R3[[r]] <- data_trt3
  
  final_theta_ctrl_R3[[r]] <- theta_ctrl
  final_theta_trt1_R3[[r]] <- theta_trt1
  final_theta_trt2_R3[[r]] <- theta_trt2
  final_theta_trt3_R3[[r]] <- theta_trt3
  
    
  end <- proc.time()
  print(end - start)
  
}  

prop_best_R1 <- prop_best_R2 <- prop_best_R3 <- rep(NA,R)
success_prob_R1 <- success_prob_R2 <- success_prob_R3 <- rep(NA,R)

clin_diff <-  0.0647669 ## omega=0.3

for(r in 1:R){
 data_all <- list(final_data_ctrl_R1[[r]],final_data_trt1_R1[[r]],
                     final_data_trt2_R1[[r]],final_data_trt3_R1[[r]])
 theta_ctrl <- final_theta_ctrl_R1[[r]]
 theta_trt1 <- final_theta_trt1_R1[[r]]
 theta_trt2 <- final_theta_trt2_R1[[r]]
 theta_trt3 <- final_theta_trt3_R1[[r]]
 
 theta_mat <- cbind(theta_trt1,theta_trt2,theta_trt3)
 theta_mean <- c(mean(theta_trt1),mean(theta_trt2),mean(theta_trt3))
 best_index <- which.max(theta_mean)
 data_best_arm <- data_all[[best_index+1]]

 prop_best_R1[r] <- prop_best_arm(data_all,data_best_arm)
 success_prob_R1[r] <- success_prob(theta_best=theta_mat[,best_index],theta_ctrl,clin_diff)

 
 data_all <- list(final_data_ctrl_R2[[r]],final_data_trt1_R2[[r]],
                  final_data_trt2_R2[[r]],final_data_trt3_R2[[r]])
 theta_ctrl <- final_theta_ctrl_R2[[r]]
 theta_trt1 <- final_theta_trt1_R2[[r]]
 theta_trt2 <- final_theta_trt2_R2[[r]]
 theta_trt3 <- final_theta_trt3_R2[[r]]
 
 theta_mat <- cbind(theta_trt1,theta_trt2,theta_trt3)
 theta_mean <- c(mean(theta_trt1),mean(theta_trt2),mean(theta_trt3))
 best_index <- which.max(theta_mean)
 data_best_arm <- data_all[[best_index+1]]
 
 prop_best_R2[r] <- prop_best_arm(data_all,data_best_arm)
 success_prob_R2[r] <- success_prob(theta_best=theta_mat[,best_index],theta_ctrl,clin_diff)
 
 data_all <- list(final_data_ctrl_R3[[r]],final_data_trt1_R3[[r]],
                  final_data_trt2_R3[[r]],final_data_trt3_R3[[r]])
 theta_ctrl <- final_theta_ctrl_R3[[r]]
 theta_trt1 <- final_theta_trt1_R3[[r]]
 theta_trt2 <- final_theta_trt2_R3[[r]]
 theta_trt3 <- final_theta_trt3_R3[[r]]
 
 theta_mat <- cbind(theta_trt1,theta_trt2,theta_trt3)
 theta_mean <- c(mean(theta_trt1),mean(theta_trt2),mean(theta_trt3))
 best_index <- which.max(theta_mean)
 data_best_arm <- data_all[[best_index+1]]
 
 prop_best_R3[r] <- prop_best_arm(data_all,data_best_arm)
 success_prob_R3[r] <- success_prob(theta_best=theta_mat[,best_index],theta_ctrl,clin_diff)
 
}

mean(prop_best_R1)
mean(prop_best_R2)
mean(prop_best_R3)

mean(success_prob_R1)
mean(success_prob_R2)
mean(success_prob_R3)

