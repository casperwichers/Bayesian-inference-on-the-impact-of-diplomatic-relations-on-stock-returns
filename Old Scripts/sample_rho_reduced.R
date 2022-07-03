sample_rho_reduced <- function(rho, alpha_star, vec_beta_star, Sigma_eta, 
                       delta, y, F_t,
                       m_min_network, m_plus_network,
                       v_min_network, v_plus_network,
                       rho_params){
  
  # Compute dimensions:
  n = nrow(y); p = ncol(y)
  # Loop over the j=1:p
  for(j in 1:p){
    # Using Beta distribution:
    if(!is.null(rho_params)){
      # Check to make sure the prior params make sense
      if(length(rho_params) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')
      
      u = (rho[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])
      
      ##################  SLICE SAMPLER WHEN USING BETA PRIOR ##################
      u = uni.slice(x0 = u, g = function(x){
        
        T = nrow(y); num_countries = ncol(y)
        rho_star = rho
        rho_star[j] = 2*x - 1
        
        A = diag(num_countries) - (diag(rho_star) %*% (delta[1]*m_min_network 
                                                       + delta[2]*m_plus_network
                                                       + delta[3]*v_min_network 
                                                       + delta[4]*v_plus_network))
        Sigma = solve(A) %*% Sigma_eta %*% t(solve(A))

        summation_term = 0
        for(t in 1:T){
          #F_tilde = x_t %*% beta_tilde
          x_temp = kronecker(t(F_t[t,]), diag(num_countries))
          mu = alpha_star + x_temp %*% vec_beta_star
          summation_term = summation_term + dmvnorm(x = y[t,], mean = mu, 
                                                    sigma = Sigma, log = TRUE)
        }
        log_llik = summation_term
        log_prior = dbeta(x, shape1 = rho_params[1], shape2 = rho_params[2], log = TRUE)
        log_dens = log_llik + log_prior
        #cat("log llik: \n", log_llik, "\n")
        #cat("log prior: \n", log_prior, "\n")
        #print(log_dens)
        return(log_dens)
        #cat("log density: \n", log_dens, "\n")
      }, w=1, m=Inf, lower = 0, upper = 1, gx0 = NULL)[1]
      rho[j] = 2*u - 1
    }
  }
  return(rho)
}
