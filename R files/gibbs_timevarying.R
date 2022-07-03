# Gibbs algorithm
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
library(ggplot2)

setwd("/Users/casper/Desktop/R files")
load("/Users/casper/Desktop/R files/6countries_192.RData")
source("uni.slice.R")
source("sample_rho_struc_tv.R")
source("log_target_dens_delta_struc_tv.R")
source("log_proposal_delta.R")
source("metropolis_delta_struc_tv.R")
source("sim_data.R")
source("plots_alphabeta.R")
source("plots_delta.R")
source("plots_everything_2layers.R")
source("plots_alphabetadelta.R")
source("plots_rho.R")
source("plots_alpha.R")
source("plots_Sigma.R")
source("norm_rows_temp_avg.R")
source("norm_rows.R")

################################################################################
################################################################################
############################  GIBBS SAMPLER ####################################
################################################################################
################################################################################
##################  SIMULATE DATA ##################
alpha_DGP = mvrnorm(n = 1, mu = rep(0,num_countries),
                    Sigma = diag(num_countries))
beta_tilde_DGP = as.vector(t(mvrnorm(n = 1, mu = rep(0,3*num_countries), 
                                     Sigma = diag(3*num_countries))))
delta_DGP = c(0.7,0.3) 
rho_DGP =  c(0.4,-0.3,0.5,0.4,-0.3,0.45)

F_t = matrix(runif(T*3), ncol=3)

Sigma_DGP = diag(runif(num_countries,0,30), num_countries)  #diag(sig2_DGP, num_countries) 

#min_network = sim_data(min_list,T)
#plus_network = sim_data(plus_list,T)
min_network = norm_rows_temp_avg(min_list)
plus_network = norm_rows_temp_avg(plus_list)

######## CHECK BONACCOLTO CONSTRAINTS #######

# # Assumption 1
# if(all(min_network[[1]] == 0) 
#    || all(plus_network[[1]] == 0)) stop("a matrix is empty")
# # Assumption 2
# if(identical(min_network[[1]], plus_network[[1]])) 
# { 
#   stop("identical matrices used")
# } 
# # Extra Assumption
# if(is.complex(eigen(min_network[[1]])$values) 
#    || is.complex(eigen(plus_network[[1]])$values))
# {
#   stop("complex eigenvalues")
# }

# ####### SAVE MIN- AND MAX EIGENVALUES FOR ALL NETWORKS ########
# m_min_evs = m_plus_evs = v_min_evs = v_plus_evs = rep(0,2)
# m_min_evs[1] = min(eigen(m_min_network[[1]])$values)
# m_min_evs[2] = max(eigen(m_min_network[[1]])$values)
# m_plus_evs[1] = min(eigen(m_plus_network[[1]])$values)
# m_plus_evs[2] = max(eigen(m_plus_network[[1]])$values)
# v_min_evs[1] = min(eigen(v_min_network[[1]])$values)
# v_min_evs[2] = max(eigen(v_min_network[[1]])$values)
# v_plus_evs[1] = min(eigen(v_plus_network[[1]])$values)
# v_plus_evs[2] = max(eigen(v_plus_network[[1]])$values)


########  CONSTRUCT DGP AND SIM RETURNS ###########
x = list()
returns = matrix(0, T, num_countries)
A_DGP = list()

for(t in 1:T){ # Since mu is time-varying (contains F_t), calc mu for every t
  A_DGP[[t]] = (diag(num_countries) - (diag(rho_DGP) %*% (delta_DGP[1]*min_network[[t]] 
                                                        + delta_DGP[2]*plus_network[[t]])))
  x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  mu_temp = solve(A_DGP[[t]]) %*% alpha_DGP + solve(A_DGP[[t]]) %*% x[[t]] %*% beta_tilde_DGP
  Sigma_temp = solve(A_DGP[[t]]) %*% Sigma_DGP %*% t(solve(A_DGP[[t]]))
  returns[t,] = mvrnorm(1, mu = mu_temp, Sigma = Sigma_temp)
}

##################  RUN THE ALGORITHM ##################
gibbs <- function(niter, burn_in, F_t, nu, rho_params, returns,
                  min_network, plus_network){
  num_countries = ncol(returns)
  T = nrow(returns)
  y = returns
  
  ################## CREATE VARIABLES ##################
  alpha = matrix(0, nrow = niter, ncol = num_countries)
  mu_alpha = rep(0,num_countries)
  Sigma_alpha = 4*diag(num_countries)
  
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  mu_beta_tilde = rep(0,3*num_countries)
  Sigma_beta_tilde = 4*diag(3*num_countries)  
  
  Sigma_eta = matrix(0, nrow = niter, ncol = num_countries)
  Sigma_eta[1,] = 1
  
  delta = matrix(0, nrow = niter, ncol = 2)
  rho = matrix(0, nrow = niter, ncol = num_countries)
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  
  ##################  DRAW INITIAL VALUES ##################
  alpha[1,] = mvrnorm(1, mu = mu_alpha, Sigma = Sigma_alpha)
  
  beta_tilde[,1] = as.vector(t(mvrnorm(n = 1, mu = rep(0,3*num_countries), 
                                       Sigma = diag(3*num_countries))))
  
  delta[1,] = rep(0.5,2)
  rho[1,] = rep(0.5, num_countries)
  
  # Setup for the algorithm, based on initial values:
  x = list()
  A = list()
  
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    A[[t]] = diag(num_countries) - (diag(rho[1,]) %*% (delta[1,1]*min_network[[t]] 
                                                     + delta[1,2]*plus_network[[t]]))
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    y_tilde[t,] = A[[t]] %*% y[t,] - x[[t]] %*% beta_tilde[,1]
  }
  
  accepts = 0
  ##################  ITERATE ##################
  for(i in 2:niter){
    print(i)
    
    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha)
    A1_alpha = T*solve(diag(Sigma_eta[i-1,]))
    An_alpha = A0_alpha + A1_alpha
    b0_alpha = solve(Sigma_alpha) %*% mu_alpha
    b1_alpha = T*solve(diag(Sigma_eta[i-1,])) %*% colMeans(y_tilde)
    bn_alpha = b0_alpha + b1_alpha
    alpha[i,] <- (mvrnorm(n=1, mu = solve(An_alpha) %*% bn_alpha,
                          Sigma = solve(An_alpha)))
    #print(alpha[i,])
    
    ######### TESTS FOR ALPHA ###########
    #alpha_star[i,] <- (mvrnorm(n=1,
    #                           mu = solve(An_alpha) %*% bn_alpha,
    #                           Sigma = solve(An_alpha)))
    #alpha[i,] = A %*% alpha_star[i,]
    #alpha_star[i,] = alpha_star_DGP
    #alpha[i,] = A_DGP %*% alpha_star[i,]
    #alpha_star[i,] = alpha_star[1,]
    #####################################
    
    # Update y_star (used in \beta_tilde) accordingly with newest \alpha
    for(t in 1:T){
      y_star[t,] = A[[t]] %*% y[t,] - alpha[i,]
    }
    
    ##### BETA TILDE #####
    A0_beta_tilde = solve(Sigma_beta_tilde)
    b0_beta_tilde = solve(Sigma_beta_tilde) %*% mu_beta_tilde
    A1_beta_tilde = 0; b1_beta_tilde = 0
    for (t in 1:T){
      A1_beta_tilde = A1_beta_tilde + t(x[[t]]) %*% solve(diag(Sigma_eta[i-1,])) %*% x[[t]] 
      b1_beta_tilde = b1_beta_tilde + t(x[[t]]) %*% solve(diag(Sigma_eta[i-1,])) %*% y_star[t,]
    }
    An_beta_tilde = A0_beta_tilde + A1_beta_tilde
    bn_beta_tilde = b0_beta_tilde + b1_beta_tilde
    beta_tilde[,i] <- (mvrnorm(n=1,
                               mu = solve(An_beta_tilde) %*% bn_beta_tilde,
                               Sigma = solve(An_beta_tilde)))
    ####### TESTS FOR BETA #########
    #beta_star_temp = matrix(vec_beta_star[,i], nrow = num_countries)
    
    #beta_bar[[i]] = A %*% beta_star_temp
    #beta_tilde[,i] = as.vector(beta_bar[[i]])
    
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
    ######################################
    
    #S_theta = matrix(data = 0, nrow = 6, ncol = 6)
    for(t in 1:T){
      # Update y_tilde (used in alpha) with newest vec_beta_star
      y_tilde[t,] = A[[t]] %*% y[t,] - x[[t]] %*% beta_tilde[,i]
    }
    
    ######### LeSage approach using one Sigma_eta ##############
    # Sigma_eta[[i]] = matrix(0, num_countries, num_countries)
    # b_tilde = 0
    # 
    # for(t in 1:T){
    #   e_temp = A[[t]] %*% y[t,] - alpha[i,] - x[[t]] %*% beta_tilde[,i]
    #   b_tilde = b_tilde + (t(e_temp) %*% e_temp)/2
    # }
    # a_tilde = T*num_countries/2
    # sig2 = rinvgamma(1,a_tilde,b_tilde) # single value
    # Sigma_eta[[i]] = diag(sig2, num_countries)
    
    ######### LeSage approach using K Sigma_eta ##############
    e_temp = matrix(0, T, num_countries)
    for(t in 1:T){
      e_temp[t,] = A[[t]] %*% y[t,] - alpha[i,] - x[[t]] %*% beta_tilde[,i]
    }
    for(j in 1:num_countries){
      b_tilde = 0
      for(t in 1:T){
        b_tilde = b_tilde + (t(e_temp[t,j]) %*% e_temp[t,j])/2
      }
      a_tilde = T/2
      Sigma_eta[i,j] = rinvgamma(1,a_tilde,b_tilde) # single value
    }
    
    ##### DELTA #####
    delta[i,] = metropolis_delta_struc_tv(start = delta[i-1,],
                                       alpha[i,], beta_tilde[,i], rho[i-1,], diag(Sigma_eta[i,]),
                                       nu, x, y,
                                       min_network, plus_network)
    #delta[i,] = delta_DGP
    #print(delta[i,])
    if(!setequal(delta[i,],delta[i-1,])){
      accepts = accepts + 1
    }
    
    ##### RHO ##### 
    #Network A will be constructed within 'sample_rho()', so no  need to update A already here
    rho[i,] = sample_rho_struc_tv(rho[i-1,], alpha[i,], beta_tilde[,i], diag(Sigma_eta[i,]),
                                delta[i,], y, x,
                                min_network, plus_network,
                                rho_params)
    #rho[i,] = rho_DGP
    #print(rho[i,])
    
    ##### UPDATE MATRIX A #####
    # Since, delta and rho are connected to A, now update A with newest values
    for(t in 1:T){
      A[[t]] = diag(num_countries) - (diag(rho[i,]) %*% (delta[i,1]*min_network[[t]]
                                                       + delta[i,2]*plus_network[[t]]))  
    }
  }
  acceptance_rate = accepts/niter
  return(list(alpha, beta_tilde, delta, rho, Sigma_eta, acceptance_rate))
}
##################  CHOOSE ALGORITHM SETUP ##################
n_trials = 10; nreps = 2000; burn_in = 500
####### SPECIFY RESULTS MATRIX ##########
results_mat = matrix(data = 0, nrow = (num_countries*2), ncol = (num_countries+3))
colnames(results_mat) <- c("Country","Measure", "alpha", "beta_1", "beta_2", "beta_3", 
                           "delta_1", "delta_2", "rho")
results_mat
for(j in 1:length(country_list)){
  results_mat[2*(j-1)+1,1] = paste0(gsub("\"", "", paste(country_list[j])))
  results_mat[2*(j-1)+2,1] = ""
  results_mat[2*(j-1)+1,2] = "Estimate"
  results_mat[2*(j-1)+2,2] = "St. Dev."
  
}

####### CHOOSE HYPERPARAMETERS #######
nu = rep(2,2)
rho_params = rep(0.5,2)

###### RUN THE ALGORITHM ########
start_time <- Sys.time()
for(i in 1:n_trials){
  list[alpha_results, beta_results, delta_results, rho_results, Sigma_results, acceptance_rate] = 
    gibbs(niter = nreps, burn_in, F_t,
          nu, rho_params, returns,
          min_network, plus_network)
  if(burn_in > 0){
    alpha_results_fil = tail(alpha_results, -burn_in)
    beta_results_fil = tail(t(beta_results), -burn_in)
    delta_results_fil = tail(delta_results, -burn_in)
    rho_results_fil = tail(rho_results, -burn_in)
    Sigma_results_fil = tail(Sigma_results, -burn_in)
  } else if(burn_in == 0){
    alpha_results_fil = alpha_results
    beta_results_fil = t(beta_results)
    delta_results_fil = delta_results
    rho_results_fil = rho_results
    Sigma_results_fil = Sigma_results
  }
  plots_everything_2layers(alpha_results_fil,beta_results_fil,
                            delta_results_fil,rho_results_fil, Sigma_results_fil)
  
  for(k in 1:num_countries){
    results_mat[2*(k-1)+1,3] = round(mean(alpha_results_fil[,k]), digits = 2) # ALPHA
    results_mat[2*(k-1)+2,3] = round(var(alpha_results_fil[,k]), digits = 2)
    results_mat[2*(k-1)+1,4] = round(mean(beta_results_fil[,k]), digits = 2) # BETA_1
    results_mat[2*(k-1)+2,4] = round(var(beta_results_fil[,k]), digits = 2)
    results_mat[2*(k-1)+1,5] = round(mean(beta_results_fil[,(k+num_countries)]), digits =2) # BETA_2
    results_mat[2*(k-1)+2,5] = round(var(beta_results_fil[,(k+num_countries)]), digits = 2)
    results_mat[2*(k-1)+1,6] = round(mean(beta_results_fil[,(k+2*num_countries)]), digits = 2) # BETA_3
    results_mat[2*(k-1)+2,6] = round(var(beta_results_fil[,(k+2*num_countries)]), digits = 2)
    results_mat[2*(k-1)+1,7] = round(mean(delta_results_fil[,1]), digits = 2) # DELTA_1
    results_mat[2*(k-1)+2,7] = round(var(delta_results_fil[,1]), digits = 2)
    results_mat[2*(k-1)+1,8] = round(mean(delta_results_fil[,2]), digits = 2) # DELTA_2
    results_mat[2*(k-1)+2,8] = round(var(delta_results_fil[,2]), digits = 2)
    results_mat[2*(k-1)+1,9] = round(mean(rho_results_fil[,k]), digits = 2) # Rho
    results_mat[2*(k-1)+2,9] = round(var(rho_results_fil[,k]), digits = 2)
  }
  # results_mat[i,1] = nu[1]
  # results_mat[i,(2:3)] = colMeans(delta_results_fil)
  # results_mat[i,(4:(ncol(results_mat)-1))] = colMeans(rho_results_fil)
  # results_mat[i,ncol(results_mat)] = rho_params[1]
  #nu = nu + 1
  #rho_params = rho_params + 1.5
}

#results_mat[nrow(results_mat),] = colMeans(results_mat[(1:nrow(results_mat)-1),])
end_time <- Sys.time()
print(end_time - start_time)
beep(sound = 11)

results_mat
xtable(results_mat, type = "latex",digits=2)
