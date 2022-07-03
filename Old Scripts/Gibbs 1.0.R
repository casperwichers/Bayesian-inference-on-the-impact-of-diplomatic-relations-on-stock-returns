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
niter = 1000
T = nrow(stock_returns)

################################################################################
################################################################################
# RUN THE ALGORITHM:
gibbs <- function(niter, burn_in, 
                  num_countries, rho = 0.5, network,
                  T = nrow(stock_returns)){ #Simple version, matrices A "known"
  ##############################################################################  
  # Setup for the algorithm:
  A = diag(num_countries) - rho*network # Assuming 1 network and 1 layer only
  y = mvrnorm(n = T, mu = rep(0,num_countries), Sigma = 4*diag(num_countries))
  F_t = matrix(runif(T*3, min =- 1, max = 1), ncol=3)
  gamma = list()
  x = list()
  for(t in 1:T){ # Construct x_t and \bar(\gamma) based on simulated F_t's and A
    x[[t]] = kronecker(t(F_t[t,]), diag(num_countries))
    gamma[[t]] = solve(A) %*% x[[t]]
  }
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
  
  y_tilde = matrix(nrow = T, ncol = num_countries)
  y_star = matrix(nrow = T, ncol = num_countries)
  ##############################################################################  
  # Draw initial values:
  alpha[1,] = mvrnorm(1, mu = mu_alpha, Sigma = Sigma_alpha) 
  beta_tilde[,1] = t(mvrnorm(1, mu=rep(0,(3*num_countries)), Sigma = 4*diag(3*num_countries)))
  Sigma[[1]] = rinvwishart(nu = nu_0, S = S_0)
  
  # Based on initial value for alpha and beta_tilde, assign values to y_tilde and y_star
  for(t in 1:T){
    y_tilde[t,] = y[t,] - solve(A) %*% x[[t]] %*% beta_tilde[,1]
    y_star[t,] = y[t,] - alpha[1,]
  }
  ##############################################################################
  # Iterate:
  for(i in 2:niter){
    ##### ALPHA #####
    A0_alpha = solve(Sigma_alpha)  #WAS WRONG FIRST
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
    A1_beta_tilde = 0
    b1_beta_tilde = 0
    for (t in 1:T){
      A1_beta_tilde = A1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% gamma[[t]]
      b1_beta_tilde = b1_beta_tilde + t(gamma[[t]]) %*% solve(Sigma[[i-1]]) %*% y_star[t,]
    }
    b0_beta_tilde = solve(Sigma_beta_tilde) %*% mu_beta_tilde
    An_beta_tilde = A0_beta_tilde + A1_beta_tilde
    bn_beta_tilde = b0_beta_tilde + b1_beta_tilde
    
    # #A1_beta_tilde = T^2 * t(gamma_bar) %*% solve(Sigma[[i-1]]) %*% gamma_bar # CHECK FOR T^2; ASSUMPTION BECAUSE 2 TIMES GAMMA BAR INSTEAD OF SUM
    # A1_beta_tilde = 
    # An_beta_tilde = A0_beta_tilde + A1_beta_tilde
    # b0_beta_tilde = solve(Sigma_beta_tilde) %*% mu_beta_tilde
    # bn_beta_tilde = b0_beta_tilde + b1_beta_tilde

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
      mu_Sigma[,t] = alpha[i,] + solve(A) %*% x[[t]] %*% beta_tilde[,i]
      S_theta[[i]] = S_theta[[i]] + (y[t,] - mu_Sigma[,t]) %*% t(y[t,] - mu_Sigma[,t]) 
      #print(S_theta[[i]])
    }
    # Add small correction value:
    S_theta[[i]] = S_theta[[i]] + diag(10^(-5), num_countries)

    ##### SIGMA #####
    S_n[[i]] = S_0 + S_theta[[i]]
    Sigma[[i]] =  rinvwishart(nu = nu_0 + T,
                              S = solve(S_n[[i]]))
  }
  return(list(colMeans(tail(alpha,-burn_in)),
                       rowMeans(tail(beta_tilde,-burn_in)), Sigma))
}

# Run the algorithm:
# list[alpha_results, beta_tilde_results, Sigma_results] = 
#   gibbs(niter = 1000, num_countries, rho = 0.5, network = v_min_list[[10]],
#       T = nrow(stock_returns))

n_trials = 20
alphas = matrix(0, nrow=n_trials,ncol=num_countries)
for(i in 1:30){
  list[temp_alpha_results, temp_beta_tilde_results, temp_Sigma_results] = 
    gibbs(niter = 1000, num_countries, rho = 0.5, network = v_min_list[[10]],
          T = nrow(stock_returns))
  alphas[i,] = tail(temp_alpha_results,1)
}

par(mfrow=c(2,2))
plot(alphas[,1], main = "Plot for alpha_1")
plot(alphas[,2], main = "Plot for alpha_2")
plot(alphas[,3], main = "Plot for alpha_3")
plot(alphas[,4], main = "Plot for alpha_4")
plot(alphas[,5], main = "Plot for alpha_5")
plot(alphas[,6], main = "Plot for alpha_6")




