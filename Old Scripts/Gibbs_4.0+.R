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
alpha_DGP = mvrnorm(n = 1, mu = c(0,0,0,0,0,0), Sigma = diag(num_countries))
beta_bar_DGP = t(mvrnorm(n = 3, mu = rep(0,num_countries), Sigma = diag(num_countries)))
delta_DGP = c(0.4,0.3,0.2,0.1) #Arbitrary, can choose and see if algo converges to this value
rho_DGP =  c(0.2,0.3,0.5,0.4,0.3,0.3)

F_t = matrix(runif(T*3), ncol=3)
Sigma_eta = diag(num_countries)#*5

m_min_network = sim_data(m_min_rows,1)
m_plus_network = sim_data(m_plus_rows,1)
v_min_network = sim_data(v_min_rows,1)
v_plus_network = sim_data(v_plus_rows,1)

####### SAVE MIN- AND MAX EIGENVALUES FOR ALL NETWORKS ########
# m_min_evs = m_plus_evs = v_min_evs = v_plus_evs = rep(0,2)
# m_min_evs[1] = min(eigen(m_min_network[[1]])$values)
# m_min_evs[2] = max(eigen(m_min_network[[1]])$values)
# m_plus_evs[1] = min(eigen(m_plus_network[[1]])$values)
# m_plus_evs[2] = max(eigen(m_plus_network[[1]])$values)
# v_min_evs[1] = min(eigen(v_min_network[[1]])$values)
# v_min_evs[2] = max(eigen(v_min_network[[1]])$values)
# v_plus_evs[1] = min(eigen(v_plus_network[[1]])$values)
# v_plus_evs[2] = max(eigen(v_plus_network[[1]])$values)
###############################################################

x = list(); #Sigma = list() # Create variables
returns = matrix(0, T, num_countries)

A_DGP = (diag(num_countries) - (diag(rho_DGP) %*% 
                                (delta_DGP[1]*m_min_network[[1]] 
                               + delta_DGP[2]*m_plus_network[[1]]
                               + delta_DGP[3]*v_min_network[[1]] 
                               + delta_DGP[4]*v_plus_network[[1]])))
Sigma_DGP = solve(A_DGP) %*% Sigma_eta %*% t(solve(A_DGP))
#print(Sigma_DGP)

beta_tilde_DGP = as.vector(beta_bar_DGP)
vec_beta_star_DGP = as.vector(solve(A_DGP) %*% beta_bar_DGP)

for(t in 1:T){ # Since mu is time-varying (contains F_t), calc mu for every t
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu_temp = alpha_DGP + x[[t]] %*% vec_beta_star_DGP # And use this mu to sample returns
  returns[t,] = mvrnorm(1, mu = mu_temp, Sigma = Sigma_DGP)
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
  
  vec_beta_star = matrix(nrow = 3*num_countries, ncol = niter)
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  
  mu_vec_beta_star = rep(0,3*num_countries)
  Sigma_vec_beta_star = 4*diag(3*num_countries)  
  
  Sigma = list(); #S_theta = list(); #S_n = list()
  nu_0 = num_countries+2  # originally: num_countries^2+10
  S_0 = 1.5*diag(num_countries) #play with, originally 1.5
  
  delta = matrix(0, nrow = niter, ncol = 4)
  rho = matrix(0, nrow = niter, ncol = num_countries)
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  
  ##################  DRAW INITIAL VALUES ##################
  alpha[1,] = mvrnorm(1, mu = mu_alpha, Sigma = Sigma_alpha) 
  Sigma[[1]] = 1.5*diag(num_countries)
  delta[1,] = c(0.1,0.05,0.05,0.8) #rep(0.25,4)
  rho[1,] = rep(0.5, num_countries)
  
  # Setup for the algorithm, based on initial values:
  x = list(); #gamma = list()
  A = diag(num_countries) - (diag(rho[1,]) %*% (delta[1,1]*m_min_network 
                                              + delta[1,2]*m_plus_network
                                              + delta[1,3]*v_min_network 
                                              + delta[1,4]*v_plus_network))
  
  beta_bar_ini = t(mvrnorm(n = 3, mu = rep(0,num_countries), Sigma = diag(num_countries)))
  vec_beta_star[,1] = as.vector(solve(A) %*% beta_bar_ini)
  beta_tilde[,1] = as.vector(beta_bar_ini)
  
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    y_tilde[t,] = y[t,] - x[[t]] %*% vec_beta_star[,1]
  }
  
  ##################  ITERATE ##################
  for(i in 2:niter){
    #print(i)
    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha)
    A1_alpha = T*solve(Sigma[[i-1]])
    An_alpha = A0_alpha + A1_alpha
    b0_alpha = solve(Sigma_alpha) %*% mu_alpha
    b1_alpha = T*solve(Sigma[[i-1]]) %*% colMeans(y_tilde)
    bn_alpha = b0_alpha + b1_alpha
    alpha[i,] <- (mvrnorm(n=1,
                          mu = solve(An_alpha) %*% bn_alpha,
                          Sigma = solve(An_alpha)))
    #alpha[i,] = alpha_DGP
    #cat("Alpha \n", alpha[i,], "\n")
    
    # Update y_star (used in \beta_tilde) accordingly with newest \alpha
    for(k in 1:num_countries){
      y_star[,k] = y[,k] - alpha[i,k]
    }
    
    ##### BETA TILDE #####
    A0_vec_beta_star = solve(Sigma_vec_beta_star)
    b0_vec_beta_star = solve(Sigma_vec_beta_star) %*% mu_vec_beta_star
    A1_vec_beta_star = 0; b1_vec_beta_star = 0
    for (t in 1:T){
      A1_vec_beta_star = A1_vec_beta_star + t(x[[t]]) %*% solve(Sigma[[i-1]]) %*% x[[t]]
      b1_vec_beta_star = b1_vec_beta_star + t(x[[t]]) %*% solve(Sigma[[i-1]]) %*% y_star[t,]
    }
    An_vec_beta_star = A0_vec_beta_star + A1_vec_beta_star
    bn_vec_beta_star = b0_vec_beta_star + b1_vec_beta_star
    vec_beta_star[,i] <- (mvrnorm(n=1,
                               mu = solve(An_vec_beta_star) %*% bn_vec_beta_star,
                               Sigma = solve(An_vec_beta_star)))
    #vec_beta_star[,i] = vec_beta_star_DGP
    beta_star_temp = matrix(vec_beta_star[,i], nrow = num_countries)
    beta_bar_temp = A %*% beta_star_temp
    #beta_bar_temp = A_DGP %*% beta_star_temp
    beta_tilde[,i] = as.vector(beta_bar_temp)

    ####### tests for beta #########
    #cat("vec_beta_star[,i] \n", vec_beta_star[,i], "\n")
    #print("beta_star_temp:")
    #print(beta_star_temp)
    #cat("beta_bar_temp \n", beta_bar_temp, "\n")
    #print(as.vector(beta_bar_temp))
    #beta_tilde[,i] = beta_tilde_ini
    #cat("Beta \n", beta_tilde[,i], "\n")
    
    #beta_star_temp = matrix(vec_beta_star[,i], nrow = num_countries)
    #beta_bar_temp = A_DGP %*% beta_star_temp
    #beta_tilde[,i] = as.vector(beta_bar_temp)
    
    #vec_beta_star[,i] = vec_beta_star_DGP
    ###############################
    
    S_theta = matrix(data = 0, nrow = 6, ncol = 6)
    for(t in 1:T){
      # Update y_tilde (used in alpha) with newest vec_beta_star
      y_tilde[t,] = y[t,] - x[[t]] %*% vec_beta_star[,i]
      ##### SIGMA #####
      mu_t_temp = alpha[i,] + x[[t]] %*% vec_beta_star[,i]
      S_theta = S_theta + (y[t,] - mu_t_temp) %*% t(y[t,] - mu_t_temp)
    }
    #print(S_theta)
    S_theta = S_theta + diag(10^(-5), num_countries) # Add small correction value
    S_n = S_0 + S_theta
    Sigma[[i]] = solve( rwish(v = nu_0 + T, S = solve(S_n)) ) # Exact same syntax as Hoff
    #print(Sigma[[i]])
    #Sigma_eta = A %*% Sigma[[i]] %*% t(A)
    #print(Sigma_eta)
    #Sigma[[i]] = Sigma_DGP
    
    ##### DELTA #####
    delta[i,] = metropolis_delta(start = delta[i-1,],
                                 beta_tilde[,i], rho[i-1,], Sigma[[i]],
                                 nu, F_t, y_star,
                                 m_min_network, m_plus_network,
                                 v_min_network, v_plus_network)
    #delta[i,] = delta_DGP
    #print(paste("DELTA:", delta[i,]))
    
    ##### RHO ##### 
    # Network A will be constructed within 'sample_rho()', so no need to update A already here
    rho[i,] = sample_rho(rho[i-1,], alpha[i,], beta_tilde[,i], Sigma[[i]],
                        delta[i,], y_star, F_t,
                        m_min_network, m_plus_network,
                        v_min_network, v_plus_network,
                        rho_params)
    #print(paste("RHO", rho[i]))
    #rho[i,] = rho_DGP
    #rho[i,] = rep(0.2,num_countries)
    
    ##### UPDATE MATRIX A #####
    # Since, delta and rho are connected to A, now update A with newest values
    A = diag(num_countries) - (diag(rho[i,]) %*% (delta[i,1]*m_min_network
                                                + delta[i,2]*m_plus_network
                                                + delta[i,3]*v_min_network
                                                + delta[i,4]*v_plus_network))
    
    # for(t in 1:T){ # and also update all terms that contain A directly
    #    y_tilde[t,] = y[t,] - x[[t]] %*% vec_beta_star[,i]
    #    y_tilde[t,] = y[t,] - (solve(A) %*% x[[t]] %*% beta_tilde[,i])
    # }
    
  }
  return(list(alpha, vec_beta_star, delta, rho))
}
##################  CHOOSE ALGORITHM SETUP ##################
n_trials = 1; nreps = 1000; burn_in = 100
nu = rep(0.5,4)
start_time <- Sys.time()
for(i in 1:n_trials){
  list[alpha_results, beta_results, delta_results, rho_results] = 
    gibbs(niter = nreps, burn_in, F_t,
          nu, rho_params = c(2,2), returns,
          m_min_network[[1]], m_plus_network[[1]],
          v_min_network[[1]], v_plus_network[[1]])
  if(burn_in > 0){
    alpha_results_fil = tail(alpha_results, -burn_in)
    beta_results_fil = tail(t(beta_results), -burn_in)
    delta_results_fil = tail(delta_results, -burn_in)
    rho_results_fil = tail(rho_results, -burn_in)
  } else if(burn_in == 0){
    alpha_results_fil = alpha_results
    beta_results_fil = t(beta_results)
    delta_results_fil = delta_results
    rho_results_fil = rho_results
  }
  end_time <- Sys.time()
  beep(sound = 11)
}
  
##################  PLOTS FOR EVERYTHING ##################
plots_everything <- function(){
  par(mfrow=c(6,3))
  for(i in 1:3){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
  for(i in 1:3){ # DELTA PLOTS
    plot(delta_results_fil[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results_fil[,i])-0.1,
                max(delta_results_fil[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
  }
  for(i in 1:3){ # DELTA PLOTS
    hist(delta_results_fil[,i], freq = F,
         breaks = seq(min(delta_results_fil[,i]) - 0.1,
                      max(delta_results_fil[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i)), "nu = ", nu[1]),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
    #acf(delta_results_fil[,i], lag.max = 10, type = "correlation", plot = TRUE)
  }
  for(i in 1:6){ # RHO PLOTS
    hist(rho_results_fil[,i], freq = F, col = "blue",
         breaks = seq(min(rho_results_fil[,i]) - 0.1, 
                      max(rho_results_fil[,i]) + 0.1, by=0.005),
         main = paste0("Histogram of results of rho", i))
    abline(v = rho_DGP[i], col="red")
  }
  for(i in 1:3){ # BETA PLOTS
    hist(Trim(beta_results_fil[,i], trim = 0.01), freq = F, col = "blue",
         breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.01) - 3),
                      max(Trim(beta_results_fil[,i],trim = 0.01) + 3), by=0.1),
         main = paste("Hist of sampled values of beta tilde", i))
    abline(v = vec_beta_star_DGP[i], col = "red")
  }
  mtext(paste("n_iter =",nreps, ", burn_in =", burn_in,
              ", start value = (0.25,0.25,0.25,0.25)",
              ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
        side = 3, line = -26, outer = TRUE)
}
################  PLOTS FOR ALPHA & BETA ##################
plots_alphabeta <- function(){
  par(mfrow=c(6,3))
  for(i in 1:6){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
  for(i in 1:12){ # BETA PLOTS
    hist(Trim(beta_results_fil[,i], trim = 0.02), freq = F, col = "blue",
         # breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.02) - 0.2),
         #              max(Trim(beta_results_fil[,i],trim = 0.02) + 0.2), by=0.1),
         breaks = 20,
         main = paste("Hist of sampled values of vec(beta_star)", i))
    abline(v = vec_beta_star_DGP[i], col = "red")
  }
}
############# PLOTS FOR ALPHA, BETA & DELTA ###############
plots_alphabetadelta <- function(){
  par(mfrow=c(6,3))
  for(i in 1:5){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
  for(i in 1:5){ # BETA PLOTS
    hist(Trim(beta_results_fil[,i], trim = 0.05), freq = F, col = "blue",
         breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.05) - 0.2),
                      max(Trim(beta_results_fil[,i],trim = 0.05) + 0.2), by=0.1),
         main = paste("Hist of sampled values of vec(beta_star)", i))
    abline(v = vec_beta_star_DGP[i], col = "red")
  }
  for(i in 1:4){
    plot(delta_results[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results[,i])-0.1,
                max(delta_results[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
    hist(delta_results_fil[,i], freq = F,
         breaks = seq(min(delta_results_fil[,i]-0.05),
                      max(delta_results_fil[,i])+0.05, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
  }
}
# ############################# PLOTS FOR ALPHA ##################################
plots_alpha <- function(){
  par(mfrow=c(6,2))
  for(i in 1:ncol(alpha_results)){
    plot(alpha_results_fil[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
         ylab = "Average value",
         col = "blue",
         ylim=c(min(alpha_results[,i])-0.5,
                max(alpha_results[,i])+0.5))
    abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
}
# ############################# PLOTS FOR DELTA ##################################
plots_delta <- function(){
  par(mfrow=c(4,2))
  for(i in 1:4){
    plot(delta_results[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results[,i])-0.1,
                max(delta_results[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
    hist(delta_results[,i], freq = F,
         breaks = seq(min(delta_results[,i]) - 0.1,
                      max(delta_results[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
  }
  mtext(paste("n_iter =",nreps, ", burn_in =", burn_in,
              ", start value = (0.25,0.25,0.25,0.25)",
              ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
        side = 3, line = -18, outer = TRUE)
}
# ############################# PLOTS FOR RHO ##################################
plots_rho <- function(){
  par(mfrow=c(6,2))
  for(i in 1:num_countries){
    plot(rho_results_fil[,i], type = "l")
    abline(h = rho_DGP[i], col = "red")
    hist(rho_results_fil[,i], freq = F,
         breaks = seq(min(rho_results_fil[,i]) - 0.1,
                      max(rho_results_fil[,i]) + 0.1, by=0.02),
         main = paste0("Histogram of results of rho", i))
    abline(v = rho_DGP[i], col="red")
  } 
}


print(end_time - start_time)                                                      




# plots_alpha()
# plots_alphabeta()
# plots_alphabetadelta()
# plots_delta()
# plots_rho()
plots_everything()




