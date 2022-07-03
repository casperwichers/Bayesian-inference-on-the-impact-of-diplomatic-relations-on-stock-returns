plots_everything_realdata <- function(alpha, beta, delta, rho, Sigma){
  par(mfrow=c(6,4))
  for(i in 1:4){ # ALPHA PLOTS
    hist(Trim(alpha[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha[,i], trim = 0.05)-0.1),
                      max(Trim(alpha[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("alpha_", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "black")
    if(exists("alpha_DGP")){
      abline(v = alpha_DGP[i], col="red")  
    }
  }
  for(i in 1:4){ # BETA PLOTS
    hist(Trim(beta[,i], trim = 0.01), freq = F, col = "black",
         breaks = seq(min(Trim(beta[,i], trim = 0.01) - 3),
                      max(Trim(beta[,i],trim = 0.01) + 3), by=0.1),
         main = paste("Hist of sampled values of beta_tilde", i))
    if(exists("beta_tilde_DGP")){
      abline(v = beta_tilde_DGP[i], col = "red")
    }
  }
  for(i in 1:2){ # DELTA TRACE PLOTS
    plot(delta[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "black",
         ylim=c(min(delta[,i])-0.1,
                max(delta[,i])+0.1))
    if(exists("delta_DGP")){
      abline(h = delta_DGP[i], col="red")
    }
  }
  for(i in 1:2){ # DELTA HIST PLOTS
    hist(delta[,i], freq = F,
         breaks = seq(min(delta[,i]) - 0.1,
                      max(delta[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i)), "nu = ", nu[1]),
         xlab = "Value",
         col = "black")
    if(exists("delta_DGP")){
      abline(v = delta_DGP[i], col="red")
    }
    #acf(delta_results_fil[,i], lag.max = 10, type = "correlation", plot = TRUE)
  }
  for(i in 1:6){ # RHO PLOTS
    hist(rho[,i], freq = F, col = "black",
         breaks = seq(min(rho[,i]) - 0.1, 
                      max(rho[,i]) + 0.1, by=0.005),
         main = paste0("Histogram of results of rho", i))
    if(exists("rho_DGP")){
      abline(v = rho_DGP[i], col="red")
    }
  }
  for(i in 1:6){
    # hist(Trim(Sigma[,i], trim = 0.05), freq = F,
    #      breaks = seq(min(Trim(Sigma[,i], trim = 0.05)-0.1),
    #                   max(Trim(Sigma[,i], trim = 0.05)+0.1), by=0.05),
    #      main = paste("Hist of MH values for",gsub(" ", "", paste("σ^2_", i)),"trimmed by 0.01"),
    hist(Sigma[,i], freq = F,
         breaks = seq(min(Sigma[,i]-10),
                      max(Sigma[,i]+10), by=0.5),
         main = paste("Hist of MH values for",gsub(" ", "", paste("σ^2_", i)),"trimmed by 0.01"),
         #xlim = c(0,80),
         xlab = "Value",
         col = "black")
    if(exists("Sigma_DGP")){
      abline(v = Sigma_DGP[i,i], col="red")
    }
  }
  # mtext(paste("# Iterations =", nreps, ", # Burned-in samples =", burn_in,
  #             "\n start value = (0.5,0.5)",
  #             ", alpha =",gsub(" ", "", paste("(",nu[1],",",nu[2],")")),
  #             "\n delta means = ", 
  #             gsub(" ", "", paste("(",round(mean(delta[,1]),digits = 3),
  #                                 ", ",round(mean(delta[,2]),digits = 3),")"))),
  #       side = 1, adj = 1, at = 3)
}