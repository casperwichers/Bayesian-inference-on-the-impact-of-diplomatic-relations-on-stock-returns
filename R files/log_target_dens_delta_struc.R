log_target_dens_delta_struc <- function(delta, alpha, beta_tilde, Sigma_eta,
                                            rho, F_t, nu, y,
                                            min_network, plus_network){
  T = nrow(y); num_countries = ncol(y)
  
  A = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network 
                                          + delta[2]*plus_network))
  # log_summation_term = 0
  # for(t in 1:T){
  #   x_temp = kronecker(t(F_t[t,]), diag(num_countries))
  #   e_temp = A %*% y[t,] - alpha - x_temp %*% beta_tilde
  #   log_summation_term = log_summation_term + ((t(e_temp) %*% e_temp))
  # }

  summation_term = 0
  for(t in 1:T){
    x_temp = kronecker(t(F_t[t,]), diag(num_countries))
    e_temp = A %*% y[t,] - alpha - x_temp %*% beta_tilde
    summation_term = summation_term + (t(e_temp) %*% solve(Sigma_eta) %*% e_temp)
  }
  llik = det(A)^T %*% summation_term^(-num_countries/2)
   
  # summation_term = 0
  # for(t in 1:T){
  #   x_temp = kronecker(t(F_t[t,]), diag(num_countries))
  #   e_temp = A %*% y[t,] - alpha - x_temp %*% beta_tilde
  #   summation_term = summation_term + (t(e_temp) %*% e_temp)/Sigma_eta
  # }

  #llik = det(A)^T %*% summation_term
  ##############
  
  #print(T*log(det(A)))
  #print((num_countries/2)*log_summation_term)
  
  #log_llik = T*log(det(A)) - summation_term/2
  log_prior = ddirichlet_R(x = delta, alpha = nu, log = FALSE)
  #print(log_llik)
  #print(log_prior)
  log
  log_posterior = log_prior + llik
  return(log_posterior)
}


