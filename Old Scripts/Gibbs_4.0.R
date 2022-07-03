# MH Algorithm with gibbs steps for delta and included, 
# Multiple rho's
# time-invariant networks 
# Update A after every iteration with newest rho and delta
##################  IMPORT PACKAGES ##################
cat("\014") 
library(MASS)
library(DirichletReg)
library(LaplacesDemon)
library(dplyr)
library(tidyr)
library(faux)
library(wordspace)
library(boot)
library(matlib)
library(gsubfn)
library(DescTools)
library(beepr)
library(MCMCpack)
##################  READ IN DATA ##################
stock_returns = matrix(unlist(read.csv("/Users/casper/Desktop/Scriptie/csv files/returns.csv")),
                       ncol = 6)
num_countries = ncol(stock_returns)
T = nrow(stock_returns)
##################  SIMULATE DATA ##################
alpha_ini = mvrnorm(n = 1, mu = c(0,0,0,0,0,0), Sigma = diag(num_countries))
beta_tilde_ini = as.vector(mvrnorm(n = 3, mu = rep(0,num_countries), Sigma = diag(num_countries)))
delta_ini = c(0.4,0.3,0.02,0.1) #Arbitrary, can choose and see if algo converges to this value
rho_ini =  c(0.2,0.3,0.5,0.7,0.6,0.8)

F_t = matrix(runif(T*3), ncol=3)
Sigma_eta = diag(num_countries)#*5

m_min_network = sim_data(m_min_rows,1)
m_plus_network = sim_data(m_plus_rows,1)
v_min_network = sim_data(v_min_rows,1)
v_plus_network = sim_data(v_plus_rows,1)

A = list(); mu = list();x = list(); Sigma = list() # Create variables
returns = matrix(0, T, num_countries)
for(t in 1:T){ # Since mu time-varying (contains F_t), calc "mu" for all t
  A[[t]] = (diag(num_countries) - (diag(rho_ini) %*% (delta_ini[1]*m_min_network[[1]] 
                                                    + delta_ini[2]*m_plus_network[[1]]
                                                    + delta_ini[3]*v_min_network[[1]] 
                                                    + delta_ini[4]*v_plus_network[[1]])))
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu[[t]] = alpha_ini + solve(A[[t]]) %*% x[[t]] %*% beta_tilde_ini # And use this mu to sample returns
  Sigma[[t]] = solve(A[[t]]) %*% Sigma_eta %*% t(solve(A[[t]]))
  Sigma_ini = Sigma[[t]]
  returns[t,] = mvrnorm(1, mu = mu[[t]], Sigma = Sigma[[t]])
}
##################  RUN THE ALGORITHM ##################
gibbs <- function(niter, burn_in, F_t, nu, rho_params, returns,
                  m_min_network, m_plus_network,
                  v_min_network, v_plus_network){
  num_countries = ncol(returns)
  T = nrow(returns)
  y = returns
  ################## CREATE VARIABLES ##################
  alpha = matrix(0,nrow = niter,ncol = num_countries)
  mu_alpha = rep(0,num_countries)
  Sigma_alpha = 4*diag(num_countries)
  
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  mu_beta_tilde = rep(0,3*num_countries)
  Sigma_beta_tilde = 4*diag(3*num_countries)
  
  Sigma = list(); S_theta = list(); S_n = list()
  mu_Sigma = matrix(nrow = niter, ncol = num_countries)
  nu_0 = num_countries+2  # originally: num_countries^2+10
  S_0 = 1.5*diag(num_countries) #play with, originally 1.5
  
  delta = matrix(0, nrow = niter, ncol = 4)
  rho = matrix(0, nrow = niter, ncol = num_countries)
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  ##################  DRAW INITIAL VALUES ##################
  alpha[1,] = mvrnorm(1, mu = mu_alpha, Sigma = Sigma_alpha) 
  beta_tilde[,1] = t(mvrnorm(1, mu=rep(0,(3*num_countries)), 
                             Sigma = Sigma_beta_tilde))
  #Sigma[[1]] = rinvwishart(nu = nu_0, S = S_0)
  Sigma[[1]] = 1.5*diag(num_countries)
  delta[1,] = rep(0.25,4)
  rho[1,] = rep(0.5, num_countries)
  
  # Setup for the algorithm, based on initial values:
  gamma = list(); x = list()
  A = diag(num_countries) - (diag(rho[1,]) %*% (delta[1,1]*m_min_network 
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
  ##################  ITERATE ##################
  for(i in 2:niter){
    #print(i)
    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha)
    A1_alpha = T*solve(Sigma[[i-1]]) # This kind of explodes
    An_alpha = A0_alpha + A1_alpha
    b0_alpha = solve(Sigma_alpha) %*% mu_alpha
    b1_alpha = T*solve(Sigma[[i-1]]) %*% colMeans(y_tilde)
    bn_alpha = b0_alpha + b1_alpha
    alpha[i,] <- (mvrnorm(n=1,
                          mu = solve(An_alpha) %*% bn_alpha,
                          Sigma = solve(An_alpha)))
    #alpha[i,] = alpha_ini
    #cat("Alpha \n", alpha[i,], "\n")
    
    # Update y_star (used in \beta_tilde) accordingly with newest \alpha
    for(t in 1:T){
      y_star[t,] = y[t,] - alpha[i,]
    }
    ##### BETA TILDE #####
    # A0_beta_tilde = solve(Sigma_beta_tilde)
    # b0_beta_tilde = solve(Sigma_beta_tilde) %*% mu_beta_tilde
    # A1_beta_tilde = 0; b1_beta_tilde = 0
    # for (t in 1:T){
    #   A1_beta_tilde = A1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% gamma[[t]]
    #   b1_beta_tilde = b1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% y_star[t,]
    # }
    # An_beta_tilde = A0_beta_tilde + A1_beta_tilde
    # bn_beta_tilde = b0_beta_tilde + b1_beta_tilde
    # beta_tilde[,i] <- (mvrnorm(n=1,
    #                            mu = solve(An_beta_tilde) %*% bn_beta_tilde,
    #                            Sigma = solve(An_beta_tilde)))
    beta_tilde[,i] = beta_tilde_ini
    #cat("Beta \n", beta_tilde[,i], "\n")
    
    S_theta[[i]] = matrix(data = 0, nrow = 6, ncol = 6)
    for(t in 1:T){
      # Change y_tilde (used in \alpha) accordingly with newest \beta_tilde
      y_tilde[t,] = y[t,] - (solve(A) %*% x[[t]] %*% beta_tilde[,i])
      # And also change mu_Sigma using the newest values of \alpha and \beta_tilde
      # next 2 lines: S_theta = sum_{t=1}^T (y_t - mu_Sigma_t) t(y_t - mu_Sigma_t)
      mu_Sigma[i,] = alpha[i,] + solve(A) %*% x[[t]] %*% beta_tilde[,i]
      S_theta[[i]] = S_theta[[i]] + (y[t,] - mu_Sigma[i,]) %*% t(y[t,] - mu_Sigma[i,])
    }
    #print(S_theta[[i]])
    # Add small correction value:
    S_theta[[i]] = S_theta[[i]] + diag(10^(-5), num_countries)
    ##### SIGMA #####
    S_n[[i]] = S_0 + S_theta[[i]]
    #print(solve(S_theta[[i]]))
    Sigma[[i]] = solve( rwish(v = nu_0 + T, S = solve(S_n[[i]])) ) # Exact same syntax as Hoff
    #### tests for Sigma ####
    # print(Sigma[[i]])
    # print(rInvWishart(df = nu_0 + T, scale = S_n[[i]]))
    # print(rinvwishart(nu = nu_0 + T, S = S_n[[i]]))
    # Sigma[[i]] =  rinvwishart(nu = nu_0 + T,
    #                           S = S_n[[i]])
    # #cat("Sigma: \n", Sigma[[i]])
    #Sigma[[i]] = Sigma_ini
    
    
    
    ##### DELTA #####
    # nreps=2, since first rep is initial value -> only one "draw" is made here
    delta[i,] = metropolis_delta(start = delta[i-1,],
                                 alpha[i,], beta_tilde[,i], rho[i-1,], Sigma[[i]],
                                 nu, F_t, returns,
                                 m_min_network, m_plus_network,
                                 v_min_network, v_plus_network)
    #delta[i,] = delta_ini
    #print(paste("DELTA:", delta[i,]))
    
    ##### RHO ##### 
    # Network A will be constructed within 'sample_rho()', so no need to update A already here
    # rho[i,] = sample_rho(rho[i-1,], alpha[i,], beta_tilde[,i], Sigma[[i]],
    #                     delta[i,], returns, F_t,
    #                     m_min_network, m_plus_network,
    #                     v_min_network, v_plus_network,
    #                     rho_params)
    #print(paste("RHO", rho[i]))
    rho[i,] = rho_ini
    
    ##### UPDATE A MATRIX #####
    # Since, delta and rho are connected to A, now update A with newest values
    A = diag(num_countries) - (diag(rho[i,]) %*% (delta[i,1]*m_min_network 
                                                + delta[i,2]*m_plus_network
                                                + delta[i,3]*v_min_network 
                                                + delta[i,4]*v_plus_network))
    for(t in 1:T){ # and also update all terms that contain A directly
      gamma[[t]] = solve(A) %*% x[[t]]
      y_tilde[t,] = y[t,] - solve(A) %*% x[[t]] %*% beta_tilde[,i]
    }
  }
  return(list(alpha, beta_tilde, delta, rho))
}
##################  CHOOSE ALGORITHM SETUP ##################
n_trials = 1; nreps = 2500; burn_in = 1000
nu = rep(1,4)
start_time <- Sys.time()
for(i in 1:n_trials){
  list[alpha_results, beta_tilde_results, delta_results, rho_results] = 
    gibbs(niter = nreps, burn_in, F_t, 
          nu, rho_params = c(1,6), returns,
          m_min_network[[1]], m_plus_network[[1]],
          v_min_network[[1]], v_plus_network[[1]])
  if(burn_in > 0){
    alpha_results_fil = tail(alpha_results, -burn_in)
    beta_tilde_results_fil = tail(t(beta_tilde_results), -burn_in)
    delta_results_fil = tail(delta_results, -burn_in)
    rho_results_fil = tail(rho_results, -burn_in)
  } else if(burn_in == 0){
    alpha_results_fil = alpha_results
    beta_tilde_results_fil = t(beta_tilde_results)
    delta_results_fil = delta_results
    rho_results_fil = rho_results
  }
  end_time <- Sys.time()
  beep(sound = 11)
  ##################  PLOTS FOR EVERYTHING ##################
  par(mfrow=c(5,3))
  for(i in 1:3){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_ini[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_ini[i], col="red")
  }
  for(i in 1:3){ # DELTA PLOTS
    plot(delta_results_fil[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results_fil[,i])-0.1,
                max(delta_results_fil[,i])+0.1))
    abline(h = delta_ini[i], col="red")
  }
  for(i in 1:3){ # DELTA PLOTS
    hist(delta_results_fil[,i], freq = F,
         breaks = seq(min(delta_results_fil[,i]) - 0.1,
                      max(delta_results_fil[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i)), "nu = ", nu[1]),
         xlab = "Value",
         col = "blue")
    abline(v = delta_ini[i], col="red")
    #acf(delta_results_fil[,i], lag.max = 10, type = "correlation", plot = TRUE)
  }
  for(i in 1:3){ # RHO PLOTS
    hist(rho_results_fil[,i], freq = F, col = "blue",
         breaks = seq(min(rho_results_fil[,i]) - 0.1, 
                      max(rho_results_fil[,i]) + 0.1, by=0.005),
         main = paste0("Histogram of results of rho", i))
    abline(v = rho_ini[i], col="red")
  }
  for(i in 1:3){ # BETA TILDE PLOTS
    hist(Trim(beta_tilde_results_fil[,i], trim = 0.01), freq = F, col = "blue",
         breaks = seq(min(Trim(beta_tilde_results_fil[,i], trim = 0.01) - 3),
                max(Trim(beta_tilde_results_fil[,i],trim = 0.01) + 3), by=0.1),
         main = paste("Hist of sampled values of beta tilde", i))
    abline(v = beta_tilde_ini[i], col = "red")
  }
  nu = nu+0.5
}
# ############################# PLOTS FOR ALPHA ##################################
par(mfrow=c(6,2))
for(i in 1:ncol(alpha_results)){
  plot(alpha_results_fil[,i], type = "l", lwd = 3,
       main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
       ylab = "Average value",
       col = "blue",
       ylim=c(min(alpha_results[,i])-0.5,
              max(alpha_results[,i])+0.5))
  abline(h = alpha_ini[i], col="red")
  hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
       breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                    max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
       main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
       xlab = "Value",
       col = "blue")
  abline(v = alpha_ini[i], col="red")
}

# ############################# PLOTS FOR DELTA ##################################
par(mfrow=c(4,2))
for(i in 1:4){
  plot(delta_results[,i], type = "l", lwd = 3,
              main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
              ylab = "Value",
              col = "blue",
              ylim=c(min(delta_results[,i])-0.1,
                     max(delta_results[,i])+0.1))
         abline(h = delta_ini[i], col="red")
  hist(delta_results[,i], freq = F,
       breaks = seq(min(delta_results[,i]) - 0.1,
                    max(delta_results[,i]) + 0.1, by=0.01),
       main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
       xlab = "Value",
       col = "blue")
  abline(v = delta_ini[i], col="red")
}
mtext(paste("n_iter =",nreps, ", burn_in =", burn_in,
            ", start value = (0.25,0.25,0.25,0.25)",
            ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
      side = 3, line = -18, outer = TRUE)

# ############################# PLOTS FOR RHO ##################################
# par(mfrow=c(6,2))
# for(i in 1:num_countries){
#   plot(rho_results_fil[,i], type = "l")
#   hist(rho_results_fil[,i], freq = F,
#        breaks = seq(min(rho_results_fil[,i]) - 0.1,
#                   max(rho_results_fil[,i]) + 0.1, by=0.02),
#        main = paste0("Histogram of results of rho", i))
#   abline(v = rho_ini[i], col="red")
# }

print(end_time - start_time)                                                      






