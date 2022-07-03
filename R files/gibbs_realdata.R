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
library(xtable)

setwd("/Users/casper/Desktop/TVP-SAR")
load("/Users/casper/Desktop/R files/6countries_192.RData")
source("uni.slice.R")
source("sample_rho_struc_tv.R")
source("log_target_dens_delta_struc_tv.R")
source("log_proposal_delta.R")
source("metropolis_delta_struc_tv.R")
source("sim_data.R")
source("plots_alphabeta.R")
source("plots_delta.R")
source("plots_everything_realdata.R")
source("plots_alphabetadelta.R")
source("plots_rho.R")
source("plots_alpha.R")
source("plots_Sigma.R")
source("norm_rows.R")
source("max_row_norm.R")
source("norm_rows_temp_avg.R")
source("plots_delta_thesis.R")
source("plots_beta.R")
source("plots_thesis.R")
source("plots_beta.R")
source("plots_delta_bw.R")

################################################################################
################################################################################
############################  GIBBS SAMPLER ####################################
################################################################################
################################################################################

start_date = min(dates)
end_date = max(dates)

################## IF USING REAL DATA ###################
stock_returns = read.csv("/Users/casper/Desktop/Scriptie/csv files/returns.csv")
stock_returns = subset(stock_returns[(stock_returns[,1] >= start_date) &
                                       (stock_returns[,1] < end_date),])
rownames(stock_returns) = stock_returns$X; stock_returns<-stock_returns[,-1]
stock_returns = stock_mat
stock_returns = stock_returns*100 # CHECK IF TRUE

global_factors = read.csv("/Users/casper/Desktop/Scriptie/csv files/global_factors.csv")
global_factors = subset(global_factors[(global_factors[,1] >= start_date) &
                                         (global_factors[,1] < end_date),])
rownames(global_factors) = global_factors$X; global_factors <- global_factors[,-1]

num_countries = ncol(stock_returns)
T = nrow(stock_returns)

returns = stock_returns
F_t = global_factors

######## NORMALIZE THE ROWS ###############
###### NORMAL ROW-NORMALIZATION ###########
min_network = norm_rows(min_list)
plus_network = norm_rows(plus_list)
###########################################
####### MAX-ROW NORMALIZATION #############
#min_network = max_row_norm(min_list)
#plus_network = max_row_norm(plus_list)
###########################################
####### DYNAMIC TEMPORAL AVERAGES #########
#min_network = norm_rows_temp_avg(min_list)
#plus_network = norm_rows_temp_avg(plus_list)

########## IF USING SIMULATED DATA ##################
# alpha_DGP = mvrnorm(n = 1, mu = rep(0,num_countries),
#                     Sigma = diag(num_countries))
# beta_tilde_DGP = as.vector(t(mvrnorm(n = 1, mu = rep(0,3*num_countries), 
#                                      Sigma = diag(3*num_countries))))
# delta_DGP = c(0.65,0.35) 
# rho_DGP =  c(0.4,0.35,0.5,-0.4,0.35,-0.3)
# 
# F_t = matrix(runif(T*3), ncol=3)
# 
# Sigma_DGP = diag(runif(num_countries,0,5), num_countries)
# x = list()
# returns = matrix(0, T, num_countries)
# A_DGP = list()
# 
# for(t in 1:T){ # Since mu is time-varying (contains F_t), calc mu for every t
#   A_DGP[[t]] = (diag(num_countries) - (diag(rho_DGP) %*% (delta_DGP[1]*min_network[[t]] 
#                                                           + delta_DGP[2]*plus_network[[t]])))
#   x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
#   mu_temp = solve(A_DGP[[t]]) %*% alpha_DGP + solve(A_DGP[[t]]) %*% x[[t]] %*% beta_tilde_DGP
#   Sigma_temp = solve(A_DGP[[t]]) %*% Sigma_DGP %*% t(solve(A_DGP[[t]]))
#   returns[t,] = mvrnorm(1, mu = mu_temp, Sigma = Sigma_temp)
# }

######## CHECK BONACCOLTO CONSTRAINTS #######

for(t in 1:T){
  # Assumption 1
  if(all(min_network[[t]] == 0)
     || all(plus_network[[t]] == 0)) stop("a matrix is empty")
  # Assumption 2
  if(identical(min_network[[t]], plus_network[[t]]))
  {
    stop("identical matrices used")
  }
  # Extra Assumption
  if(is.complex(eigen(min_network[[t]])$values)
     || is.complex(eigen(plus_network[[t]])$values))
  {
    print("complex eigenvalues")
  }
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
  Sigma_alpha = diag(num_countries)*2#*0.1
  
  beta_tilde = matrix(nrow = 3*num_countries, ncol = niter)
  mu_beta_tilde = rep(0,3*num_countries)
  Sigma_beta_tilde = diag(3*num_countries)*2
  
  a_sigma_prior = 3#9/4
  b_sigma_prior = 2#5/4
  
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
  
  x = list()
  A = list()
  
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    A[[t]] = diag(num_countries) - (diag(rho[1,]) %*% (delta[1,1]*min_network[[t]] 
                                                     + delta[1,2]*plus_network[[t]]))
    x[[t]] = kronecker(t(as.vector(unlist(F_t[t,]))), diag(num_countries))
    y_tilde[t,] = A[[t]] %*% unlist(y[t,]) - x[[t]] %*% beta_tilde[,1]
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

    # Update y_star (used in \beta_tilde) accordingly with newest \alpha
    for(t in 1:T){
      y_star[t,] = A[[t]] %*% unlist(y[t,]) - alpha[i,]
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
    
    for(t in 1:T){
      # Update y_tilde (used in alpha) with newest vec_beta_star
      y_tilde[t,] = A[[t]] %*% unlist(y[t,]) - x[[t]] %*% beta_tilde[,i]
    }
    
    ######### LeSage approach using K Sigma_eta ##############
    e_temp = matrix(0, T, num_countries)
    for(t in 1:T){
      e_temp[t,] = A[[t]] %*% unlist(y[t,]) - alpha[i,] - x[[t]] %*% beta_tilde[,i]
    }
    for(j in 1:num_countries){
      b_tilde = 0
      for(t in 1:T){
        b_tilde = b_tilde + (t(e_temp[t,j]) %*% e_temp[t,j])/2
      }
      a_tilde = (T/2)+1
      Sigma_eta[i,j] = rinvgamma(1,a_tilde,b_tilde) # single value
      #Sigma_eta[i,j] = rinvgamma(1,a_tilde+a_sigma_prior,b_tilde+b_sigma_prior) # single value
      #Sigma_eta[i,j] = rinvgamma(1,a_tilde+3,b_tilde+2) # single value
    }
    
    ##### DELTA #####
    delta[i,] = metropolis_delta_struc_tv(start = delta[i-1,],
                                          alpha[i,], beta_tilde[,i], rho[i-1,], diag(Sigma_eta[i,]),
                                          nu, x, y,
                                          min_network, plus_network)
    if(!setequal(delta[i,],delta[i-1,])){
      accepts = accepts + 1
    }
    
    ##### RHO ##### 
    #Network A will be constructed within 'sample_rho()', so no  need to update A already here
    rho[i,] = sample_rho_struc_tv(rho[i-1,], alpha[i,], beta_tilde[,i], diag(Sigma_eta[i,]),
                                  delta[i,], y, x,
                                  min_network, plus_network,
                                  rho_params)
    print(rho[i,])
    
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
n_trials = 1; nreps = 10000; burn_in = 2500
####### SPECIFY RESULTS MATRIX ##########
results_mat = matrix(data = 0, nrow = (num_countries*3), ncol = (num_countries+2))
colnames(results_mat) <- c("Country","Measure", "$alpha$", "$beta_{oil}$", "$beta_{mat}$", "$beta_{GEPU}$", 
                           "$Sigma_eta$", #"delta_1", "delta_2", 
                           "$rho$")
country_list2 = country_list
country_list2[1] = "UK";country_list2[4] = "Russia";country_list2[6] = "US"
results_mat
for(j in 1:length(country_list)){
  results_mat[3*(j-1)+1,1] = paste0(gsub("\"", "", paste(country_list2[j])))
  results_mat[3*(j-1)+2,1] = ""
  results_mat[3*(j-1)+3,1] = ""
  results_mat[3*(j-1)+1,2] = "Post. Mean"
  results_mat[3*(j-1)+2,2] = "SD"
  results_mat[3*(j-1)+3,2] = "CI"
}
results_mat
#xtable(results_mat, type = "latex",digits=2)

####### CHOOSE HYPERPARAMETERS #######
nu = rep(2,2)
rho_params = rep(0.5,2)

start_time <- Sys.time()
for(i in 1:n_trials){
  ##### SIMULATE DATA ########
  # alpha_DGP = mvrnorm(n = 1, mu = rep(0,num_countries),
  #                     Sigma = diag(num_countries))
  # beta_tilde_DGP = as.vector(t(mvrnorm(n = 1, mu = rep(0,3*num_countries),
  #                                      Sigma = diag(3*num_countries))))
  # delta_DGP = c(0.65,0.35)
  # rho_DGP =  runif(6,0.1,0.8)
  # 
  # F_t = matrix(runif(T*3,-1,1), ncol=3)
  # 
  # Sigma_DGP = diag(runif(num_countries,0,5), num_countries)
  # x = list()
  # returns = matrix(0, T, num_countries)
  # A_DGP = list()
  # 
  # for(t in 1:T){ # Since mu is time-varying (contains F_t), calc mu for every t
  #   A_DGP[[t]] = (diag(num_countries) - (diag(rho_DGP) %*% (delta_DGP[1]*min_network[[t]]
  #                                                           + delta_DGP[2]*plus_network[[t]])))
  #   x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
  #   mu_temp = solve(A_DGP[[t]]) %*% alpha_DGP + solve(A_DGP[[t]]) %*% x[[t]] %*% beta_tilde_DGP
  #   Sigma_temp = solve(A_DGP[[t]]) %*% Sigma_DGP %*% t(solve(A_DGP[[t]]))
  #   returns[t,] = mvrnorm(1, mu = mu_temp, Sigma = Sigma_temp)
  # }
  ###### RUN THE ALGORITHM ########
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
  #plots_everything_realdata(alpha_results_fil,beta_results_fil,
  #                      delta_results_fil,rho_results_fil, Sigma_results_fil)
  plots_rho(rho_results_fil, country_list2)
  plots_beta(beta_results_fil,country_list2)
  plots_country(alpha_results_fil, beta_results_fil, delta_results_fil,
                 rho_results_fil, sigma_results_fil, country_list2)
  # plots_delta_bw(delta_results_fil)
  plots_delta_thesis(delta_results_fil)
  
  cr = 1.96
  
  for(k in 1:num_countries){
    results_mat[3*(k-1)+1,3] = round(mean(alpha_results_fil[,k]), digits = 2) # ALPHA
    results_mat[3*(k-1)+2,3] = gsub(" ", "",paste("(",round(sd(alpha_results_fil[,k]), digits = 2),")"))
    results_mat[3*(k-1)+3,3] = gsub(" ", "",paste("[",round(quantile(alpha_results_fil[,k], p=0.05/2),digits=2),",",
                                                  round(quantile(alpha_results_fil[,k],p=1-0.05/2),digits=2),"]"))
    results_mat[3*(k-1)+1,4] = round(mean(beta_results_fil[,k]), digits = 2) # BETA_1
    results_mat[3*(k-1)+2,4] = gsub(" ", "",paste("(",round(sd(beta_results_fil[,k]), digits = 2),")"))
    results_mat[3*(k-1)+3,4] = gsub(" ", "",paste("[",round(quantile(beta_results_fil[,k], p=0.05/2),digits=2),",",
                                                  round(quantile(beta_results_fil[,k],p=1-0.05/2),digits=2),"]"))
    results_mat[3*(k-1)+1,5] = round(mean(beta_results_fil[,(k+num_countries)]), digits =2) # BETA_2
    results_mat[3*(k-1)+2,5] = gsub(" ", "",paste("(",round(sd(beta_results_fil[,(k+num_countries)]), digits = 2),")"))
    results_mat[3*(k-1)+3,5] = gsub(" ", "",paste("[",round(quantile(beta_results_fil[,(k+num_countries)], p=0.05/2),digits=2),",",
                                                  round(quantile(beta_results_fil[,(k+num_countries)],p=1-0.05/2),digits=2),"]"))
    results_mat[3*(k-1)+1,6] = round(mean(beta_results_fil[,(k+2*num_countries)]), digits = 2) # BETA_3
    results_mat[3*(k-1)+2,6] = gsub(" ", "",paste("(",round(sd(beta_results_fil[,(k+2*num_countries)]), digits = 2),")"))
    results_mat[3*(k-1)+3,6] = gsub(" ", "",paste("[",round(quantile(beta_results_fil[,(k+2*num_countries)], p=0.05/2),digits=2),",",
                                                  round(quantile(beta_results_fil[,(k+2*num_countries)],p=1-0.05/2),digits=2),"]"))
    results_mat[3*(k-1)+1,7] = round(mean(Sigma_results_fil[,k]), digits = 2) # Sigma
    results_mat[3*(k-1)+2,7] = gsub(" ", "",paste("(",round(sd(Sigma_results_fil[,k]), digits = 2),")"))
    results_mat[3*(k-1)+3,7] = gsub(" ", "",paste("[",round(quantile(Sigma_results_fil[,k], p=0.05/2),digits=2),",",
                                                  round(quantile(Sigma_results_fil[,k],p=1-0.05/2),digits=2),"]"))
    results_mat[3*(k-1)+1,8] = round(mean(rho_results_fil[,k]), digits = 2) # Rho
    results_mat[3*(k-1)+2,8] = gsub(" ", "",paste("(",round(sd(rho_results_fil[,k]), digits = 2),")"))
    results_mat[3*(k-1)+3,8] = gsub(" ", "",paste("[",round(quantile(rho_results_fil[,k], p=0.05/2),digits=2),",",
                                                  round(quantile(rho_results_fil[,k],p=1-0.05/2),digits=2),"]"))
    
  }
  assign(paste0("results_mat",i), results_mat)
}

end_time <- Sys.time()
print(end_time - start_time)
beep(sound = 11)

results_mat

xtable(results_mat, type = "latex",digits=2)

quantile(delta_results_fil[,1], p = 0.05/2)
quantile(delta_results_fil[,1], p = 1-0.05/2)
quantile(delta_results_fil[,2], p = 0.05/2)
quantile(delta_results_fil[,2], p = 1-0.05/2)


