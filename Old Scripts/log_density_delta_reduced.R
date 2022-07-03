log_target_density_delta_reduced <- function(delta, alpha, beta_tilde, Sigma_eta,
                                     rho, F_t, nu, y,
                                     m_min_network, m_plus_network,
                                     v_min_network, v_plus_network){
  T = nrow(y); num_countries = ncol(y)
  
  A = diag(num_countries) - diag(rho) %*% (delta[1]*m_min_network 
                                           + delta[2]*m_plus_network
                                           + delta[3]*v_min_network 
                                           + delta[4]*v_plus_network)
  # CHANGES ARE MADE HERE
  Sigma = solve(A) %*% Sigma_eta %*% t(solve(A))

  alpha_star = alpha
  vec_beta_star = beta_tilde
  
  summation_term = 0
  for(t in 1:T){
    #F_tilde = x_t %*% beta_tilde
    x_temp = kronecker(t(F_t[t,]), diag(num_countries))
    #mu_temp = y[t,] - (alpha_star + x_temp %*% vec_beta_star)
    #summation_term = (summation_term + (t(mu_temp) %*% solve(Sigma) %*% (mu_temp)))
    mu = alpha_star + x_temp %*% vec_beta_star
    summation_term = summation_term + dmvnorm(x = y[t,], mean = mu, 
                                              sigma = Sigma, log = TRUE)
  }
  log_llik = summation_term
  #log_llik = -(T/2)*log(det(Sigma)) - summation_term/2
  
  log_prior = ddirichlet_R(x = delta, alpha = nu, log = TRUE)
  # log_prior = ((nu[1]-1)*log(delta[1]) + (nu[2]-1)*log(delta[2])
  #            + (nu[3]-1)*log(delta[3]) + (nu[4]-1)*log(delta[4]))
  
  log_posterior = log_prior + log_llik
  return(log_posterior)
}