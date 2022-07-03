metropolis_doublestep <- function(start, rho,
                             alpha, beta_tilde, 
                             Sigma_eta, y, F_t, 
                             m_min_network, m_plus_network,
                             v_min_network, v_plus_network,
                             nu, rho_params){
  delta_old = matrix(data = start, nrow = 1, ncol = 4)
  rho_old = rho

  repeat {
    theta_star <- rdirichlet(1, alpha = nu)
    if(max(theta_star) < 0.8 && min(theta_star) > 0.05) break
  }
  
  rho_star = sample_rho(rho_old, alpha, beta_tilde,
                        Sigma_eta, delta_old, y, F_t,
                        m_min_network, m_plus_network,
                        v_min_network, v_plus_network,
                        rho_params)

  # print("alpha:")
  # print(alpha)
  # print("alpha_DGP:")
  # print(alpha_DGP)
  # print("beta_tilde:")
  # print(beta_tilde)
  # print("beta tilde DGP:")
  # print(beta_tilde_DGP)
  
  #alpha = alpha_DGP
  #beta_tilde = beta_tilde_DGP
  #theta_star = matrix(delta_DGP,1,4)
  #rho_star = rho_DGP
  # print("Hier komen de samples:")
  # print(theta_star)
  # print(rho_star)
  # print("dit zijn de oude samples:")
  # print(delta_old)
  # print(rho_old)
  
  ############# COMPUTE A AND B ################
  
  b = (log_target_density_delta(delta = theta_star, alpha, beta_tilde, Sigma_eta,
                                rho = rho_star, F_t, nu, y,
                                m_min_network, m_plus_network,
                                v_min_network, v_plus_network) 
       + log_proposal_delta(delta_old, nu))
  print(log_target_density_delta(delta = theta_star, alpha, beta_tilde, Sigma_eta,
                                 rho = rho_star, F_t, nu, y,
                                 m_min_network, m_plus_network,
                                 v_min_network, v_plus_network))
  print(log_proposal_delta(delta_old, nu))
  cat("new density", b, "\n")
  a = (log_target_density_delta(delta = delta_old, alpha, beta_tilde, Sigma_eta,
                                rho = rho_old, F_t, nu, y,
                                m_min_network, m_plus_network,
                                v_min_network, v_plus_network)
       + log_proposal_delta(theta_star, nu))
  print(log_target_density_delta(delta = delta_old, alpha, beta_tilde, Sigma_eta,
                                  rho = rho_old, F_t, nu, y,
                                  m_min_network, m_plus_network,
                                  v_min_network, v_plus_network))
  print(log_proposal_delta(theta_star, nu))
  
  cat("old density", a, "\n")
  H = min(1,exp(b-a)) # Take ratio like this since using logs
  #cat("H:", round(H), "\n")
  runif = runif(1)
  if(runif < H){ 
    theta <- theta_star
    new_rho = rho_star
    cat("Accepted! \n old theta: \n", delta_old, "\n", "new theta: \n", theta_star, "\n")
    cat("The old posterior density was:", a, "\n")
    cat("And the new posterior density is:", b, "\n")
    cat("H = ", H, "\n")
    cat("The uniform draw was: \n", runif, "\n")
    pause(10)
  } else {
    theta <- delta_old
    new_rho <- rho_old
  }
  return(list(theta,new_rho))
}