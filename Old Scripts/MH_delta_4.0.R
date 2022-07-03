# MH Algorithm for delta 
# time-invariant networks

library(MASS)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(faux)
library(wordspace)
library(boot)
library(DirichletReg)
library(matlib)
library(gsubfn)
library(beepr)
cat("\014") 

##################  TARGET DENSITY ##################
log_target_density_delta <- function(delta, alpha, beta_tilde, Sigma,
                                     rho, returns, F_t, nu, 
                                     m_min_network, m_plus_network,
                                     v_min_network, v_plus_network){
  T = nrow(returns); num_countries = ncol(returns)
  A = list() #;Sigma = list();
  x = list(); F_tilde = list();
  
  y_tilde = matrix(nrow = T, ncol = num_countries) #y_tilde = y_t - alpha
  for(i in 1:ncol(returns)){
    y_tilde[,i] = returns[,i]-alpha[i]
  }
  
  summation_term = 0
  for(t in 1:T){
    A[[t]] = diag(num_countries) - diag(rho) %*% (delta[1]*m_min_network 
                                                + delta[2]*m_plus_network
                                                + delta[3]*v_min_network 
                                                + delta[4]*v_plus_network)
    # Sigma[[t]] = solve(A[[t]]) %*% Sigma_eta %*% t(solve(A[[t]]))
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    F_tilde[[t]] = x[[t]] %*% beta_tilde
    temp = y_tilde[t,] - (solve(A[[t]]) %*% F_tilde[[t]])
    summation_term = (summation_term + (t(temp) %*% solve(Sigma) %*% (temp)))
  }
  log_llik = log(det(Sigma)) - summation_term/2
  
  # log_prior_term = ddirichlet(x = matrix(data=delta, ncol = 4), 
  #                             alpha = nu, log = TRUE)
  log_prior = ((nu[1]-1)*log(delta[1]) 
             + (nu[2]-1)*log(delta[2])
             + (nu[3]-1)*log(delta[3]) 
             + (nu[4]-1)*log(delta[4]))
  
  log_posterior = log_prior + log_llik
  return(log_posterior) #Single value
}

##################  PROPOSAL DENSITY ##################
log_proposal_delta <- function(value, nu){ #not smapling, evaluating
  log_dens = ddirichlet_R(x = matrix(data=value, ncol = 4),
                        alpha = nu, log = TRUE) #hyperparams always same
  manual_proposal = (nu[1]-1)*log(delta[1]) 
                + (nu[2]-1)*log(delta[2])
                + (nu[3]-1)*log(delta[3]) 
                + (nu[4]-1)*log(delta[4])
  cat("Proposal with ddirichlet:", log_dens)
  cat("Proposal manually:", manual_proposal)
  
  return(log_dens)
}

##################  MH ALGORITHM TO BE USED IN GIBBS SAMPLER ##################
metropolis_delta <- function(start, nreps = 2,burn_in,
                             alpha, beta_tilde, rho, Sigma,
                             nu, F_t, returns,
                             m_min_network, m_plus_network,
                             v_min_network, v_plus_network){
  theta <- matrix(nrow = nreps, ncol = 4)
  theta[1,] <- start
  #print(paste("theta[1,]:", theta[1,]))
  H = rep(0,nreps)
  for (i in 2:nreps){
    theta_star <- rdirichlet(1, alpha = nu) # <- this will be updated metropolis_rw_delta
    #print(theta_star)
    b = (log_target_density_delta(delta = theta_star, alpha, beta_tilde, Sigma,
                                  rho, returns, F_t, nu,
                                  m_min_network, m_plus_network,
                                  v_min_network, v_plus_network) 
         + log_proposal_delta(theta[i-1,], nu))
    a = (log_target_density_delta(delta = theta[i-1,], alpha, beta_tilde, Sigma,
                                  rho, returns, F_t, nu,
                                  m_min_network, m_plus_network,
                                  v_min_network, v_plus_network) 
         + log_proposal_delta(theta_star, nu))
    #cat("b: \n", b, "\n")
    #cat("a: \n", a, "\n")
    #cat("b-a: \n", b-a, "\n")
    #cat("exp(b-a): \n", exp(b-a), "\n")
    
    H[i] = min(1,exp(b-a)) # Take ratio like this since using logs
    #print(H[i])
    if(runif(1) < H[i]){ 
      theta[i,] <- theta_star
    } else {
      theta[i,] <- start
    }
  }
  return(tail(theta,1))
}
