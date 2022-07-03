# MH Algorithm for delta 
# time-invariant networks

library(MASS)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(faux)
library(wordspace)
library(boot)
#library(DirichletReg)
library(matlib)
library(gsubfn)
library(beepr)
cat("\014") 

##################  TARGET DENSITY ##################
log_target_density_delta <- function(delta, alpha, vec_beta_star, Sigma,
                                     rho, F_t, nu, y_star,
                                     m_min_network, m_plus_network,
                                     v_min_network, v_plus_network){
  T = nrow(y_star); num_countries = ncol(y_star)
  
  A = diag(num_countries) - diag(rho) %*% (delta[1]*m_min_network 
                                           + delta[2]*m_plus_network
                                           + delta[3]*v_min_network 
                                           + delta[4]*v_plus_network)
  summation_term = 0
  for(t in 1:T){
    #F_tilde = x_t %*% beta_tilde
    #F_tilde_temp = kronecker(t(F_t[t,]), diag(num_countries)) %*% beta_tilde
    #temp = y_star[t,] - (solve(A) %*% F_tilde_temp)
    
    temp = y_star[t,] - kronecker(t(F_t[t,]), diag(num_countries)) %*% vec_beta_star
    
    summation_term = (summation_term + (t(temp) %*% solve(Sigma) %*% (temp)))
  }
  log_llik = log(det(Sigma)) - summation_term/2
  
  # log_prior_term = ddirichlet(x = matrix(data=delta, ncol = 4), 
  #                             alpha = nu, log = TRUE)
  log_prior = ((nu[1]-1)*log(delta[1]) 
               + (nu[2]-1)*log(delta[2])
               + (nu[3]-1)*log(delta[3]) 
               + (nu[4]-1)*log(delta[4]))
  #print(log_prior)
  #print(log_llik)
  log_posterior = log_prior + log_llik
  return(-log_posterior) #Single value
}

##################  PROPOSAL DENSITY ##################
log_proposal_delta <- function(delta, nu){ #not smapling, evaluating
  # log_dens = ddirichlet_R(x = matrix(data=delta, ncol = 4),
  #                         alpha = nu, log = TRUE) #hyperparams always same
  log_dens = ddirichlet_R(x = matrix(data = delta, ncol = 4), alpha = nu, log=TRUE)
  # manual_proposal = (nu[1]-1)*log(delta[1]) 
  #                 + (nu[2]-1)*log(delta[2])
  #                 + (nu[3]-1)*log(delta[3]) 
  #                 + (nu[4]-1)*log(delta[4])
  # cat("Proposal with ddirichlet: \n", log_dens, "\n")
  # cat("Proposal manually: \n", manual_proposal, "\n")
  #cat("Proposal density:", log_dens)
  return(log_dens)
}

##################  MH ALGORITHM TO BE USED IN GIBBS SAMPLER ##################
metropolis_delta <- function(start,
                             alpha, vec_beta_star, rho, Sigma,
                             nu, F_t, y_star,
                             m_min_network, m_plus_network,
                             v_min_network, v_plus_network){
  theta_ini = matrix(data = start, nrow = 1, ncol = 4)
  theta_star <- rdirichlet(1, alpha = nu)
  #print(theta_star)
  ############
  # print(theta_ini)
  # print(nrow(theta_ini))
  # print(ncol(theta_ini))
  # cat("typeof theta_ini:", typeof(theta_ini), "\n")
  # 
  # print(theta_star)
  # print(nrow(theta_star))
  # print(ncol(theta_star))
  # cat("typeof theta_star:", typeof(theta_star), "\n")
  #theta_star = c(0.4,0.3,0.2,0.1)
  # cat("Theta_ini: \n", theta_ini, "\n")
  # cat("Theta_star: \n", theta_star, "\n")
  
  # repeat {
  #   theta_star <- rdirichlet(1, alpha = nu)
  #   #theta_star <- rdirichlet(1, alpha = theta_ini*10)
  #   #print(theta_star)
  #   if(min(theta_star > 0)) break
  # }
  #print(theta_star)
  ###########
  
  b = (log_target_density_delta(delta = theta_star, alpha, vec_beta_star, Sigma,
                                rho, F_t, nu, y_star,
                                m_min_network, m_plus_network,
                                v_min_network, v_plus_network) 
       + log_proposal_delta(theta_ini, nu))
  
  a = (log_target_density_delta(delta = theta_ini, alpha, vec_beta_star, Sigma,
                                rho, F_t, nu, y_star,
                                m_min_network, m_plus_network,
                                v_min_network, v_plus_network)
       + log_proposal_delta(theta_star, nu))
  
  H = min(1,exp(b-a)) # Take ratio like this since using logs
  #cat("H:", round(H), "\n")
  runif = runif(1)
  if(runif < H){ 
    theta <- theta_star
    cat("Accepted! \n old theta: \n", theta_ini, "\n", "new theta: \n", theta_star, "\n")
    cat("The old posterior density was:", a, "\n")
    cat("And the new posterior density is:", b, "\n")
    cat("H = ", H, "\n")
    cat("The uniform draw was: \n", runif, "\n")
    
  } else {
    theta <- theta_ini
  }
  return(theta)
}


