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

#############################################################################
stock_returns = matrix(unlist(read.csv("/Users/casper/Desktop/Scriptie/csv files/returns.csv")),
                       ncol =6, nrow =)
T = nrow(stock_returns)
num_countries = ncol(stock_returns)

#############################################################################
#############################################################################

# Simulate Data
alpha = mvrnorm(n = 1, mu = c(0,0,0,0,0,0), Sigma = diag(num_countries))
beta_tilde  = as.vector(mvrnorm(n = 3, mu = rep(0,num_countries), Sigma = diag(num_countries)))
F_t = matrix(runif(T*3), ncol=3)
rho = 0.5 # Arbitrary
Sigma_eta = diag(num_countries)
#Sigma =  rinvwishart(nu = (num_countries^2+10), S = diag(num_countries))
delta = c(0.4,0.3,0.2,0.1) #Arbitrary, can choose and see if algo converges to this value

v_plus_network = sim_data(v_plus_rows,1)
v_min_network = sim_data(v_min_rows,1)
m_plus_network = sim_data(m_plus_rows,1)
m_min_network = sim_data(m_min_rows,1)

A = list(); mu = list();x = list(); Sigma = list() # Create variables
returns = matrix(0, T, num_countries)
for(t in 1:T){ # Since mu time-varying (contains F_t), calc "mu" for all t
  A[[t]] = (diag(num_countries) - (rho*(delta[1]*m_min_network[[1]] + delta[2]*m_plus_network[[1]]
                                        + delta[3]*v_min_network[[1]] + delta[4]*v_plus_network[[1]])))
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu[[t]] = alpha + solve(A[[t]]) %*% x[[t]] %*% beta_tilde # And use this mu to sample returns
  Sigma[[t]] = solve(A[[t]]) %*% Sigma_eta %*% t(solve(A[[t]]))
  returns[t,] = mvrnorm(1, mu = mu[[t]], Sigma = Sigma[[t]])
}

#############################################################################
#############################################################################

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
    A[[t]] = diag(6) - rho * (delta[1]*m_min_network + delta[2]*m_plus_network
                               +delta[3]*v_min_network + delta[4]*v_plus_network)
    # Sigma[[t]] = solve(A[[t]]) %*% Sigma_eta %*% t(solve(A[[t]]))
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    F_tilde[[t]] = x[[t]] %*% beta_tilde
    temp = y_tilde[t,] - (solve(A[[t]]) %*% F_tilde[[t]])
    summation_term = (summation_term + (t(temp) %*% solve(Sigma) %*% (temp)))
  }
  
  log_prior_term = ((nu[1]-1)*log(delta[1]) + (nu[2]-1)*log(delta[2])
                    + (nu[3]-1)*log(delta[3]) + (nu[4]-1)*log(delta[4]))
  
  log_dens = log(det(Sigma)) - summation_term/2 + log_prior_term
  return(log_dens) #Single value
}

#############################################################################
#############################################################################

log_proposal_delta <- function(value, nu){ #not smapling, evaluating
  log_dens = ddirichlet(x = matrix(data=value, ncol = 4),
                        alpha = nu, log = TRUE) #hyperparams always same
  return(log_dens)
}

#############################################################################
#############################################################################
metropolis_delta <- function(start, nreps, burn_in = 500,
                             alpha, beta_tilde, rho, Sigma,
                             nu, F_t, returns,
                             m_min_network, m_plus_network,
                             v_min_network, v_plus_network){
  theta <- matrix(nrow = nreps, ncol = 4)
  theta[1,] <- start
  H = rep(0,nreps)
  for (i in 2:nreps){
    theta_star <- rdirichlet(1, alpha = nu)
    
    b = (log_target_density_delta(delta = theta_star, alpha, beta_tilde, Sigma,
                                  rho, returns, F_t, nu,
                                  m_min_network, m_plus_network,
                                  v_min_network, v_plus_network) 
         + log_proposal_delta(theta_star,nu))
    a = (log_target_density_delta(delta = theta[i-1,], alpha, beta_tilde, Sigma,
                                  rho, returns, F_t, nu,
                                  m_min_network, m_plus_network,
                                  v_min_network, v_plus_network) 
         + log_proposal_delta(theta[i-1,],nu))
    #cat("b: \n", b, "\n")
    #cat("a: \n", a, "\n")
    #cat("b-a: \n", b-a, "\n")
    #cat("exp(b-a): \n", exp(b-a), "\n")
    
    H[i] = min(1,exp(b-a)) # Take ratio like this since using logs
    #print(H[i])
    if(runif(1) < H[i]){ 
      theta[i,] <- theta_star
    } else {
      theta[i,] <- theta[i-1,]
    }
  }
  #return(list(theta, H[-1]))
  return(tail(theta,1))
}

# n_reps=1000; burn_in_period = 0; n_trials = 1
# nu = c(1,1,1,1)
# H_mat = matrix(0, nrow = n_reps-1, ncol = n_trials)
# theta_averages = matrix(0, nrow= n_trials, ncol=4)
# #thetas = matrix(0, nrow = n_reps-burn_in_period, ncol = n_trials)
# thetas = list()
# start_time <- Sys.time()
# for(i in 1:n_trials){
#   print(i)
#   list[thetas[[i]], H_mat[,i]] = metropolis_delta(start = rep(0.25,4),
#                                                   nreps = n_reps, burn_in = burn_in_period,
#                                                   alpha, beta_tilde, rho, Sigma_eta,
#                                                   nu, F_t, returns,
#                                                   m_min_network[[1]], m_plus_network[[1]],
#                                                   v_min_network[[1]], v_plus_network[[1]])
#   if(burn_in_period > 0){
#     theta_averages[i,] = colMeans(tail(thetas[[i]],-burn_in))
#     thetas_fil[[i]] = tail(thetas[[i]],-burn_in_period)
#   } else if(burn_in_period == 0){
#     theta_averages[i,] = colMeans(thetas[[i]])
#     thetas_fil[[i]] = thetas[[i]]
#   }
# }
# 
# end_time <- Sys.time()
# end_time - start_time
# 
# par(mfrow=c(2,2))
# for(i in 1:4){  # Plot the averages of all trials
#   plot(theta_averages[,i], 
#        main = paste("Plot of averages of delta",i),
#        ylab = "Average value",
#        col = "blue",
#        ylim=c(min(theta_averages[,i])-0.1, 
#               max(theta_averages[,i])+0.1))
#   abline(h = delta[i], col="red")
# }
# mtext(paste("Plots for delta 1-4, n_iter =", n_reps, 
#             "and burn_in =", burn_in_period), 
#       side = 3, line = -18, outer = TRUE)
# 
# ################################################################################
# 
# par(mfrow=c(2,2))
# for(i in 1:4){   # Plot Histograms of separate trials
#   hist(thetas_fil[[1]][,i], freq = F,
#        breaks = seq(min(thetas_fil[[1]][,i])-0.05, 
#                     max(thetas_fil[[1]][,i])+0.05, by=0.01),
#        col = "lightblue",
#        xlab = "Parameter Value",
#        main = paste("Hist of",gsub(" ", "", paste("δ", i)),
#                     ", true value =", delta[i]))
#   abline(v = delta[i], col="red", lwd=3, lty=2)
# }
# mtext(paste("Hists for delta,", n_reps, "iterations,", 
#             burn_in_period, "burn-in samples and hyperparam alpha =",
#             gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))), 
#       side = 3, line = -18, outer = TRUE)
# 
# accprob = sum(H_mat[,1] == 1)/length(H_mat[,1])*100
# 
# par(mfrow=c(2,2))
# for(i in 1:4){ 
#   plot(thetas[[1]][,i], type = "l", lwd = 3,
#        main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
#        ylab = "Average value",
#        col = "blue",
#        ylim=c(min(thetas[[1]][,i])-0.1, 
#               max(thetas[[1]][,i])+0.1))
#   abline(h = delta[i], col="red")
# }
# 
# mtext(paste("n_iter =",n_reps, ", burn_in =", burn_in_period,
#             ", start value = (0.25,0.25,0.25,0.25)",
#             ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
#       side = 3, line = -18, outer = TRUE)
# 
# #beep(sound = 1)
# 
# 
# # metropolis_delta(start = rep(0.25,4),
# #                  nreps = 2, burn_in = 0,
# #                  alpha, beta_tilde, rho, Sigma_eta,
# #                  nu, F_t, returns,
# #                  m_min_network[[1]], m_plus_network[[1]],
# #                  v_min_network[[1]], v_plus_network[[1]])


