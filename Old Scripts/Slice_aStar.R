cat("\014") 
# Slice Sampler

#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.

uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL){
  # Here, I added that all necessary variables are also parameterized
  # Check the validity of the arguments.
  if (!is.numeric(x0) || length(x0)!=1 #<-2nd can be removed when working with multiple rho's 
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }
  
  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1
  # Find the log density at the initial point, if not already known.
  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }
  # Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)
  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  {
    x1 <- runif(1,L,R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)
    
    if (gx1>=logy) break
    
    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }
  # Return the point sampled, with its log density attached as an attribute.
  attr(x1,"log.density") <- gx1
  return (x1)
}
##################  ALGORITHM TO BE USED IN GIBBS SAMPLER ##################
sample_rho <- function(rho, alpha, beta_tilde, Sigma, 
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
        
        alpha_star = solve(A) %*% alpha
        vec_beta_star = as.vector(solve(A) %*% matrix(beta_tilde, nrow = num_countries))
        
        summation_term = 0
        for(t in 1:T){
          x_temp = kronecker(t(F_t[t,]), diag(num_countries))
          mu_temp = y[t,] - (alpha_star + x_temp %*% vec_beta_star)
          summation_term = (summation_term + (t(mu_temp) %*% solve(Sigma) %*% (mu_temp)))
        }
        
        log_llik  = log(det(Sigma)) - summation_term/2
        log_prior = dbeta(x, shape1 = rho_params[1], shape2 = rho_params[2], log = TRUE)
        log_dens = log_llik + log_prior
        #cat("log llik: \n", log_llik, "\n")
        #cat("log prior: \n", log_prior, "\n")
        return(log_dens[1,1])
        #cat("log density: \n", log_dens, "\n")
      }, w=1, m=Inf, lower = 0, upper = 1, gx0 = NULL)[1]
      rho[j] = 2*u - 1
    }
  }
  return(rho)
}
