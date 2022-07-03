# MH Algorithm with gibbs steps for delta and included, 
# Multiple rho's
# time-invariant networks 
# Update A after every iteration with newest rho and delta

##################  IMPORT PACKAGES ##################
cat("\014")
rm(list=ls())
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
library(matrixcalc)
library(mvtnorm)
library(profvis)

setwd("/Users/casper/Desktop/R files") #Set your working directory
load("/Users/casper/Desktop/R files/matrices_plus_rows.RData")
source("uni.slice.R")
source("sample_rho_reduced.R")
source("log_target_density_delta_reduced.R")
source("log_proposal_delta.R")
source("metropolis_delta_reduced.R")
source("sim_data.R")
source("plots_alpha.R")
source("plots_alphabeta.R")
source("plots_delta.R")
source("plots_everything.R")
source("plots_alphabetadelta.R")
source("plots_rho.R")

##################  READ IN DATA ##################
stock_returns = read.csv("/Users/casper/Desktop/Scriptie/csv files/returns.csv")
rownames(stock_returns) = stock_returns$X; stock_returns<-stock_returns[,-1]

num_countries = ncol(stock_returns)
T = nrow(stock_returns)

##################  SIMULATE DATA ##################
alpha_star_DGP = mvrnorm(n = 1, mu = rep(0,num_countries), 
                         Sigma = diag(num_countries))
beta_star_DGP = matrix(t(mvrnorm(n = 1, mu = rep(0,3*num_countries), 
                                 Sigma = diag(3*num_countries))), nrow = num_countries)

vec_beta_star_DGP = as.vector(beta_star_DGP)
delta_DGP = c(0.4,0.3,0.2,0.1) #Arbitrary, can choose and see if alg converges to this value
rho_DGP =  c(0.2,0.3,0.5,0.4,0.3,0.3)

F_t = matrix(runif(T*3), ncol=3)
Sigma_eta = diag(num_countries)

m_min_network = sim_data(m_min_rows,1)
m_plus_network = sim_data(m_plus_rows,1)
v_min_network = sim_data(v_min_rows,1)
v_plus_network = sim_data(v_plus_rows,1)

######## CHECK BONACCOLTO CONSTRAINTS #######

# Assumption 1
if(all(m_min_network[[1]] == 0) 
   || all(m_plus_network[[1]] == 0) 
   || all(v_min_network[[1]] == 0) 
   || all(v_plus_network[[1]] == 0)) stop("a matrix is empty")
# Assumption 2
if(identical(m_min_network[[1]], m_plus_network[[1]]) 
   || identical(m_min_network[[1]], v_min_network[[1]])
   || identical(m_min_network[[1]], v_plus_network[[1]])
   || identical(m_plus_network[[1]], v_min_network[[1]])
   || identical(m_plus_network[[1]], v_plus_network[[1]])
   || identical(v_min_network[[1]], v_plus_network[[1]]))
  { 
    stop("identical matrices used")
  } 
# Extra Assumption
if(is.complex(eigen(m_min_network[[1]])$values) 
      || is.complex(eigen(m_plus_network[[1]])$values) 
      || is.complex(eigen(v_min_network[[1]])$values)
      || is.complex(eigen(v_plus_network[[1]])$values)) 
  {
    stop("complex eigenvalues")
  }

####### SAVE MIN- AND MAX EIGENVALUES FOR ALL NETWORKS ########
m_min_evs = m_plus_evs = v_min_evs = v_plus_evs = rep(0,2)
m_min_evs[1] = min(eigen(m_min_network[[1]])$values)
m_min_evs[2] = max(eigen(m_min_network[[1]])$values)
m_plus_evs[1] = min(eigen(m_plus_network[[1]])$values)
m_plus_evs[2] = max(eigen(m_plus_network[[1]])$values)
v_min_evs[1] = min(eigen(v_min_network[[1]])$values)
v_min_evs[2] = max(eigen(v_min_network[[1]])$values)
v_plus_evs[1] = min(eigen(v_plus_network[[1]])$values)
v_plus_evs[2] = max(eigen(v_plus_network[[1]])$values)

########  CONSTRUCT DGP AND SIM RETURNS ###########
x = list()
returns = matrix(0, T, num_countries)

A_DGP = (diag(num_countries) - (diag(rho_DGP) %*% 
                                  (delta_DGP[1]*m_min_network[[1]] 
                                   + delta_DGP[2]*m_plus_network[[1]]
                                   + delta_DGP[3]*v_min_network[[1]] 
                                   + delta_DGP[4]*v_plus_network[[1]])))
Sigma_DGP = solve(A_DGP) %*% Sigma_eta %*% t(solve(A_DGP))

alpha_DGP = A_DGP %*% alpha_star_DGP
beta_bar_DGP = A_DGP %*% beta_star_DGP
beta_tilde_DGP = as.vector(beta_bar_DGP)

for(t in 1:T){ # Since mu is time-varying (contains F_t), calc mu for every t
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu_temp = alpha_star_DGP + x[[t]] %*% vec_beta_star_DGP # And use this mu to sample returns
  returns[t,] = mvrnorm(1, mu = mu_temp, Sigma = Sigma_DGP)
}

##################  RUN THE ALGORITHM ##################
gibbs <- function(switch_1, switch_2, niter, F_t, nu, rho_params, returns,
                  m_min_network, m_plus_network,
                  v_min_network, v_plus_network){
  num_countries = ncol(returns)
  T = nrow(returns)
  y = returns
  ################## CREATE VARIABLES ##################
  alpha_star = alpha = matrix(0,nrow = niter,ncol = num_countries)
  mu_alpha_star = rep(0,num_countries)
  Sigma_alpha_star = 4*diag(num_countries)
  
  vec_beta_star = matrix(nrow = 3*num_countries, ncol = niter)
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  
  mu_vec_beta_star = rep(0,3*num_countries)
  Sigma_vec_beta_star = 4*diag(3*num_countries)  
  
  Sigma = list(); #S_theta = list(); #S_n = list()
  nu_0 = num_countries^2+10  # originally: num_countries^2+10
  S_0 = 1.5*diag(num_countries) #play with, originally 1.5
  Sigma_eta = diag(num_countries)
  
  delta = matrix(0, nrow = niter, ncol = 4)
  rho = matrix(0, nrow = niter, ncol = num_countries)
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  
  ##################  DRAW INITIAL VALUES ##################
  alpha_star[1,] = mvrnorm(1, mu = mu_alpha_star, Sigma = Sigma_alpha_star)
  
  beta_star_ini = t(mvrnorm(n = 3, mu = rep(0,num_countries), 
                            Sigma = diag(num_countries)))
  vec_beta_star[,1] = as.vector(beta_star_ini)
  
  delta[1,] = rep(0.25,4)
  rho[1,] = rep(0.5, num_countries)
  
  # Setup for the algorithm, based on initial values:
  x = list()
  A = diag(num_countries) - (diag(rho[1,]) %*% (delta[1,1]*m_min_network 
                                              + delta[1,2]*m_plus_network
                                              + delta[1,3]*v_min_network 
                                              + delta[1,4]*v_plus_network))
  
  alpha[1,] = A %*% alpha_star[1,]
  beta_bar_ini = A %*% beta_star_ini
  beta_tilde[,1] = as.vector(beta_bar_ini)
  Sigma[[1]] = solve(A) %*% Sigma_eta %*% t(solve(A))
  
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    y_tilde[t,] = y[t,] - x[[t]] %*% vec_beta_star[,1]
  }
  
  ##################  ITERATE ##################
  for(i in 2:niter){
    if(i == switch_1){
      # print("Here comes the switch:")
      # print("Alpha stats:")
      # print(alpha_star_DGP)
      # print(alpha_star[(switch_1-1),])
      # print(colMeans(alpha_star[(switch_1-5):(switch_1-1),]))
      # print(colMeans(alpha_star[(switch_1-5):(switch_1-1),]))
      # print("____________________")
      # print("Beta_stats:")
      # print(head(vec_beta_star_DGP))
      # print("Last values:")
      # print(vec_beta_star[1:6,(switch_1-1)])
      # print("Rowmeans:")
      # print(rowMeans(vec_beta_star[1:6,(switch_1-10):(switch_1-1)]))
      # print("Most recent rowmeans:")
      # print(rowMeans(vec_beta_star[1:6,(switch_1-10):(switch_1-1)]))
      plots_alphabeta(alpha_star[1:(switch_1-1),], vec_beta_star[1:6,1:(switch_1-1)])
      pause(10)
    }
    if(i == switch_2){
      #print("Here comes the second switch:")
      #print(delta_DGP)
      #print(delta[(switch_2-1),])
      plots_delta(delta[switch_1:(switch_2-1),])
      #plots_alpha(alpha_star[switch_1:(switch_2-1),])
      #print("Alpha stuff:")
      #print(colMeans(alpha_star[(switch_2-5):(switch_2-1),]))
      #print("Beta stuff:")
      #print(rowMeans(vec_beta_star[,(switch_2-5):(switch_2-1)]))
      pause(10)
    }
    print(i)

    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha_star)
    A1_alpha = T*solve(Sigma[[i-1]])
    An_alpha = A0_alpha + A1_alpha
    b0_alpha = solve(Sigma_alpha_star) %*% mu_alpha_star
    b1_alpha = T*solve(Sigma[[i-1]]) %*% colMeans(y_tilde)
    bn_alpha = b0_alpha + b1_alpha
    alpha_star[i,] <- (mvrnorm(n=1,
                               mu = solve(An_alpha) %*% bn_alpha,
                               Sigma = solve(An_alpha)))
    alpha[i,] = A %*% alpha_star[i,]
    
    # Update y_star (used in \beta_tilde) accordingly with newest \alpha_star
    for(k in 1:num_countries){
      y_star[,k] = y[,k] - alpha_star[i,k]
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

    beta_star_temp = matrix(vec_beta_star[,i], nrow = num_countries)
    beta_bar_temp = A %*% beta_star_temp
    beta_tilde[,i] = as.vector(beta_bar_temp)
    
    S_theta = matrix(data = 0, nrow = 6, ncol = 6)
    for(t in 1:T){
      # Update y_tilde (used in alpha) with newest vec_beta_star
      y_tilde[t,] = y[t,] - x[[t]] %*% vec_beta_star[,i]
      ##### SIGMA #####
      mu_t_temp = alpha_star[i,] + x[[t]] %*% vec_beta_star[,i]
      S_theta = S_theta + (y[t,] - mu_t_temp) %*% t(y[t,] - mu_t_temp)
    }
    S_theta = S_theta + diag(10^(-5), num_countries) # Add small correction value
    S_n = S_0 + S_theta
    Sigma[[i]] = solve( rwish(v = nu_0 + T, S = solve(S_n)) ) # Exact same syntax as Hoff
    
    ##### DELTA #####
    if(i >= switch_1){
      delta[i,] = metropolis_delta_reduced(start = delta[i-1,],
                                   alpha_star[i,], vec_beta_star[,i], rho[i-1,], 
                                   Sigma_eta, nu, F_t, y,
                                   m_min_network, m_plus_network,
                                   v_min_network, v_plus_network)
      print(delta[i,])
    } else if(i < switch_1){
      delta[i,] = rep(0.25,4)
    }
    
    ##### RHO ##### 
    # Network A will be constructed within 'sample_rho()', so no need to update A already here
    if(i >= switch_2 ){
      rho[i,] = sample_rho_reduced(rho[i-1,], alpha_star[i,], vec_beta_star[,i], 
                           Sigma_eta, delta[i,], y, F_t,
                           m_min_network, m_plus_network,
                           v_min_network, v_plus_network,
                           rho_params)
      print(rho[i,])
    }  else if(i < switch_2){
      rho[i,] = rep(0.3,num_countries)
    }
    
    ##### UPDATE MATRIX A #####
    # Since, delta and rho are connected to A, now update A with newest values
    A = diag(num_countries) - (diag(rho[i,]) %*% (delta[i,1]*m_min_network
                                                  + delta[i,2]*m_plus_network
                                                  + delta[i,3]*v_min_network
                                                  + delta[i,4]*v_plus_network))
    
    
    # for(t in 1:T){ # and also update all terms that contain A directly
    #   updt_beta_bar = matrix(beta_tilde[,i], nrow = num_countries)
    #   updt_beta_star = solve(A) %*% updt_beta_bar
    #   y_tilde[t,] = updt_beta_star %*% F_t[t,]
    #   # y_tilde[t,] = y[t,] - (x[[t]] %*% 
    #   #               as.vector((solve(A) %*% 
    #   #                            matrix(beta_tilde[,i], nrow = num_countries))))
    #   # print(y_tilde[t,])
    # }
  }
  return(list(alpha_star, vec_beta_star, delta[switch_1:niter,], rho[switch_2:niter,]))
  #return(list(alpha_star, vec_beta_star, delta[switch_1:(switch_2-1),], rho[switch_2:niter,]))
}
##################  CHOOSE ALGORITHM SETUP ##################
n_trials = 1
burn_in = 0
nu = rep(4,4)
start_time <- Sys.time()
for(i in 1:n_trials){
  list[alpha_results, beta_results, delta_results, rho_results] = 
    gibbs(switch_1 = 2000, switch_2 = 8000, niter = 8250,
          F_t, nu, rho_params = c(2,2), returns,
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

plots_everything()

print(end_time - start_time)                                 

# bayesian estimation of spatial autogressiver models 1997 and newe ones too
# start from him and search other for other authors "SAR"
# 

