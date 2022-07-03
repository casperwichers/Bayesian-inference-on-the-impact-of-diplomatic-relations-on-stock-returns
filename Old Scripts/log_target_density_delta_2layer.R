log_target_density_delta_2layer <- function(delta, alpha, beta_tilde, Sigma_eta,
                                     rho, F_t, nu, y,
                                     min_network, plus_network){
  T = nrow(y); num_countries = ncol(y)
  
  A = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network[[1]] 
                                          + delta[2]*plus_network[[1]]))
  alpha_star = solve(A) %*% alpha
  vec_beta_star = as.vector(solve(A) %*% matrix(beta_tilde, nrow = num_countries))
  Sigma = solve(A) %*% Sigma_eta %*% t(solve(A))
  
  summation_term = 0
  for(t in 1:T){
    x_temp = kronecker(t(F_t[t,]), diag(num_countries))
    mu_temp = alpha_star + x_temp %*% vec_beta_star
    summation_term = summation_term + dmvnorm(x = y[t,], mean = mu_temp, 
                                              sigma = Sigma, log = TRUE)
  }
  log_llik = summation_term
  log_prior = ddirichlet_R(x = delta, alpha = nu, log = TRUE)
  log_posterior = log_prior + log_llik
  return(log_posterior)
}