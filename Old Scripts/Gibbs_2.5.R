# MH Algorithm with delta and rho included, 
# Time-invariant networks
# No update of A once new deltas are in

# MH Algorithm with delta and rho included
# Time-invariant networks 
# No update of A once new deltas are in

cat("\014") 
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

stock_returns = matrix(unlist(read.csv("/Users/casper/Desktop/Scriptie/csv files/returns.csv")),
                       ncol = 6)
num_countries = ncol(stock_returns)
T = nrow(stock_returns)

#############################################################################
#############################################################################
# Simulate Data
alpha_ini = mvrnorm(n = 1, mu = c(0,0,0,0,0,0), Sigma = diag(num_countries))
beta_tilde_ini  = as.vector(mvrnorm(n = 3, mu = rep(0,num_countries), Sigma = diag(num_countries)))
delta_ini = c(0.4,0.3,0.2,0.1) #Arbitrary, can choose and see if algo converges to this value
rho_ini = 0.6

F_t = matrix(runif(T*3), ncol=3)
Sigma_eta = diag(num_countries)
#rho = 0.5 # Arbitrary

v_plus_network = sim_data(v_plus_rows,1)
v_min_network = sim_data(v_min_rows,1)
m_plus_network = sim_data(m_plus_rows,1)
m_min_network = sim_data(m_min_rows,1)

A = list(); mu = list();x = list(); Sigma = list() # Create variables
returns = matrix(0, T, num_countries)
for(t in 1:T){ # Since mu time-varying (contains F_t), calc "mu" for all t
  A[[t]] = (diag(num_countries) - (rho_ini*(delta_ini[1]*m_min_network[[1]] + delta_ini[2]*m_plus_network[[1]]
                                        + delta_ini[3]*v_min_network[[1]] + delta_ini[4]*v_plus_network[[1]])))
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu[[t]] = alpha_ini + solve(A[[t]]) %*% x[[t]] %*% beta_tilde_ini # And use this mu to sample returns
  Sigma[[t]] = solve(A[[t]]) %*% Sigma_eta %*% t(solve(A[[t]]))
  returns[t,] = mvrnorm(1, mu = mu[[t]], Sigma = Sigma[[t]])
}

################################################################################
################################################################################
# RUN THE ALGORITHM:
gibbs_2.5 <- function(niter, burn_in, F_t, nu, rho_params, returns,
                  m_min_network, m_plus_network,
                  v_min_network, v_plus_network){
  
  num_countries = ncol(returns)
  T = nrow(returns)
  y = returns
  ##############################################################################  
  # Create variables used in iterations:
  alpha = matrix(0,nrow = niter,ncol = num_countries)
  mu_alpha = rep(0,num_countries)
  Sigma_alpha = 4*diag(num_countries)
  
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  mu_beta_tilde = rep(0,3*num_countries)
  Sigma_beta_tilde = 4*diag(3*num_countries)
  
  Sigma = list(); S_theta = list(); S_n = list()
  mu_Sigma = matrix(nrow = num_countries, ncol = niter)
  nu_0 = num_countries^2+10
  S_0 = 2*diag(num_countries)
  
  delta = matrix(0, nrow = niter, ncol = 4)
  head(delta)
  rho = rep(0, niter)
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  acc_prob = 0
  ##############################################################################  
  # Draw initial values:
  alpha[1,] = mvrnorm(1, mu = mu_alpha, Sigma = Sigma_alpha) 
  beta_tilde[,1] = t(mvrnorm(1, mu=rep(0,(3*num_countries)), Sigma = 4*diag(3*num_countries)))
  Sigma[[1]] = rinvwishart(nu = nu_0, S = S_0)
  delta[1,] = rep(0.25,4)
  rho[1] = runif(1,min=-1,max=1)
  
  # Setup for the algorithm, based on initial values:
  gamma = list(); x = list()
  A = diag(num_countries) - (rho[1]*(delta[1,1]*m_min_network 
                                     + delta[1,2]*m_plus_network
                                     + delta[1,3]*v_min_network 
                                     + delta[1,4]*v_plus_network))
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    gamma[[t]] = solve(A) %*% x[[t]]
    y_tilde[t,] = y[t,] - solve(A) %*% x[[t]] %*% beta_tilde[,1]
    y_star[t,] = y[t,] - alpha[1,]
  }
  ##############################################################################
  # Iterate:
  for(i in 2:niter){
    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha)
    A1_alpha = T*solve(Sigma[[i-1]])
    An_alpha = A0_alpha + A1_alpha
    b0_alpha = solve(Sigma_alpha) %*% mu_alpha
    b1_alpha = T*(solve(Sigma[[i-1]])) %*% colMeans(y_tilde)
    bn_alpha = b0_alpha + b1_alpha
    
    alpha[i,] <- (mvrnorm(n=1,
                          mu = solve(An_alpha) %*% bn_alpha,
                          Sigma = solve(An_alpha)))
    #cat("Alpha \n", alpha[i,], "\n")
    
    # Update y_star (used in \beta_tilde) accordingly with newest \alpha
    for(t in 1:T){
      y_star[t,] = y[t,] - alpha[i,]
    }
    ##### BETA TILDE #####
    A0_beta_tilde = solve(Sigma_beta_tilde)
    A1_beta_tilde = 0; b1_beta_tilde = 0
    for (t in 1:T){
      A1_beta_tilde = A1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% gamma[[t]]
      b1_beta_tilde = b1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% y_star[t,]
    }
    b0_beta_tilde = solve(Sigma_beta_tilde) %*% mu_beta_tilde
    An_beta_tilde = A0_beta_tilde + A1_beta_tilde
    bn_beta_tilde = b0_beta_tilde + b1_beta_tilde
    
    beta_tilde[,i] <- (mvrnorm(n=1,
                               mu = solve(An_beta_tilde) %*% bn_beta_tilde,
                               Sigma = solve(An_beta_tilde)))
    #cat("Beta \n", beta_tilde[,i], "\n")
    
    S_theta[[i]] = matrix(data = 0, nrow = 6, ncol = 6)
    for(t in 1:T){ 
      # Change y_tilde (used in \alpha) accordingly with newest \beta_tilde
      y_tilde[t,] = y[t,] - (solve(A) %*% x[[t]] %*% beta_tilde[,i])
      # And also change mu_Sigma using the newest values of \alpha and \beta_tilde
      # next 2 lines: S_theta = sum_{t=1}^T (y_t - mu_Sigma_t) t(y_t - mu_Sigma_t)
      mu_Sigma[,i] = alpha[i,] + solve(A) %*% x[[t]] %*% beta_tilde[,i]
      S_theta[[i]] = S_theta[[i]] + (y[t,] - mu_Sigma[,i]) %*% t(y[t,] - mu_Sigma[,i]) 
    }
    #print(S_theta[[i]])
    # Add small correction value:
    S_theta[[i]] = S_theta[[i]] + diag(10^(-5), num_countries)
    
    ##### SIGMA #####
    S_n[[i]] = S_0 + S_theta[[i]]
    Sigma[[i]] =  rinvwishart(nu = nu_0 + T,
                              S = solve(S_n[[i]]))
    
    ##### DELTA #####
    # nreps=2, since first rep is initial value -> only one "draw" is made here
    delta[i,] = metropolis_delta(start = delta[i-1,],
                                 nreps = 2, burn_in = 0,
                                 alpha[i,], beta_tilde[,i], rho[i-1], Sigma[[i]],
                                 nu, F_t, returns,
                                 m_min_network, m_plus_network,
                                 v_min_network, v_plus_network)
    if(delta[i,] != delta[i-1,]){
      acc_prob = acc_prob + 1
    }
    
    ##### RHO ##### 
    rho[i] = sample_rho(rho[i-1], alpha[i,], beta_tilde[,i], Sigma[[i]],
                        delta[i,], returns, F_t,
                        m_min_network, m_plus_network,
                        v_min_network, v_plus_network,
                        rho_params)
    
    # UPDATE A
    # Since, delta connected to A, now update A with newest values
    A = diag(num_countries) - (rho[i]*(delta[i,1]*m_min_network
                                       + delta[i,2]*m_plus_network
                                       + delta[i,3]*v_min_network
                                       + delta[i,4]*v_plus_network))
    for(t in 1:T){ # and also update all terms that contain A directly
      gamma[[t]] = solve(A) %*% x[[t]]
      y_tilde[t,] = y[t,] - solve(A) %*% x[[t]] %*% beta_tilde[,i]
    }
    
  }
  print(paste("the acceptance probability is", acc_prob/niter))
  return(list(alpha, beta_tilde, delta, rho, acc_prob))
}

n_trials = 1; nreps = 10000; burn = 0
nu = c(0.5,0.5,0.5,0.5)
for(i in 1:n_trials){
  list[alpha_results, beta_tilde_results, delta_results, rho_results, acc_prob_val] = 
    gibbs_2.5(niter = nreps, burn, F_t, nu , rho_params = c(3,2), returns,
          m_min_network[[1]], m_plus_network[[1]],
          v_min_network[[1]], v_plus_network[[1]])
}

################################################################################
# RUN ALGORITHM WITH PLOTS FOR DELTA AND ACCEPTANCE PROBABILITIES
# n_trials = 1; nreps = 7500; burn = 0
# nu = c(1,1,1,1); acc_prob_mat = matrix(data = 0, nrow = n_trials, ncol = 2)
# for(i in 1:n_trials){
#   acc_prob_mat[i,1] = nu[1] 
#   list[alpha_results, beta_tilde_results, delta_results, rho_results, acc_prob_val] = 
#     gibbs(niter = nreps, burn, F_t, nu , rho_params = c(3,2), returns,
#           m_min_network[[1]], m_plus_network[[1]],
#           v_min_network[[1]], v_plus_network[[1]])
#   acc_prob_mat[i,2] = acc_prob_val
#   
#   par(mfrow=c(2,2))
#   for(i in 1:4){ 
#     plot(delta_results[,i], type = "l", lwd = 3,
#          main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
#          ylab = "Average value",
#          col = "blue",
#          ylim=c(min(delta_results[,i])-0.1, 
#                 max(delta_results[,i])+0.1))
#     abline(h = delta_ini[i], col="red")
#   }
#   nu = nu+0.5
# }

################################################################################
############################# PLOTS FOR ALPHA ##################################
par(mfrow=c(3,2))
for(i in 1:ncol(alpha_results)){
  plot(alpha_results[,i], type = "l", lwd = 3,
       main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
       ylab = "Average value",
       col = "blue",
       ylim=c(min(alpha_results[,i])-0.5,
              max(alpha_results[,i])+0.5))
  #abline(h = alpha_ini[i], col="red")
}

# mtext(paste("n_iter =",nreps, ", burn_in =", burn,
#             ", start value = (0.25,0.25,0.25,0.25)",
#             ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
#       side = 3, line = -18, outer = TRUE)

################################################################################
############################# PLOTS FOR DELTA ##################################

par(mfrow=c(2,2))
for(i in 1:4){ 
  plot(delta_results[,i], type = "l", lwd = 3,
       main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
       ylab = "Average value",
       col = "blue",
       ylim=c(min(delta_results[,i])-0.1, 
              max(delta_results[,i])+0.1))
  abline(h = delta_ini[i], col="red")
}

mtext(paste("n_iter =",nreps, ", burn_in =", burn,
            ", start value = (0.25,0.25,0.25,0.25)",
            ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
      side = 3, line = -18, outer = TRUE)

par(mfrow=c(2,2))
for(i in 1:4){   # Plot Histograms of separate trials
  hist(delta_results[,i], freq = F,
       breaks = seq(min(delta_results[,i])-0.05,
                    max(delta_results[,i])+0.05, by=0.01),
       col = "lightblue",
       xlab = "Parameter Value",
       main = paste("Hist of",gsub(" ", "", paste("δ", i)),
                    ", true value =", delta_ini[i]))
  abline(v = delta_ini[i], col="red", lwd=3, lty=2)
}
mtext(paste("Hists for delta,", n_reps, "iterations,",
            burn_in_period, "burn-in samples and hyperparam alpha =",
            gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
      side = 3, line = -18, outer = TRUE)


plot(rho_results, ylab = "Value",
     main = "Plot of values for rho",
     col = "blue", type = "l")
