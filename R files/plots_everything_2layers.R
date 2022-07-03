plots_everything_2layers <- function(alpha, beta, delta, rho, Sigma){
  par(mfrow=c(6,4))
  for(i in 1:4){ # ALPHA PLOTS
    hist(Trim(alpha[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha[,i], trim = 0.05)-0.1),
                      max(Trim(alpha[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("alpha_", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
  for(i in 1:4){ # BETA PLOTS
    hist(Trim(beta[,i], trim = 0.01), freq = F, col = "blue",
         breaks = seq(min(Trim(beta[,i], trim = 0.01) - 3),
                      max(Trim(beta[,i],trim = 0.01) + 3), by=0.1),
         main = paste("Hist of sampled values of beta_tilde", i))
    abline(v = beta_tilde_DGP[i], col = "red")
  }
  for(i in 1:2){ # DELTA TRACE PLOTS
    plot(delta[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta[,i])-0.1,
                max(delta[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
  }
  for(i in 1:2){ # DELTA HIST PLOTS
    hist(delta[,i], freq = F,
         breaks = seq(min(delta[,i]) - 0.1,
                      max(delta[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i)), "nu = ", nu[1]),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
    #acf(delta_results_fil[,i], lag.max = 10, type = "correlation", plot = TRUE)
  }
  for(i in 1:6){ # RHO PLOTS
    hist(rho[,i], freq = F, col = "blue",
         breaks = seq(min(rho[,i]) - 0.1, 
                      max(rho[,i]) + 0.1, by=0.005),
         main = paste0("Histogram of results of rho", i))
    abline(v = rho_DGP[i], col="red")
  }
  for(i in 1:6){
    hist(Trim(Sigma[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(Sigma[,i], trim = 0.05)-0.1),
                      max(Trim(Sigma[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("σ^2_", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = Sigma_DGP[i,i], col="red")
  }
  mtext(paste("# Iterations =", nreps, ", # Burned-in samples =", burn_in,
              "\n start value = (0.5,0.5)",
              ", alpha =",gsub(" ", "", paste("(",nu[1],",",nu[2],")")),
              "\n delta means = ", 
              gsub(" ", "", paste("(",round(mean(delta[,1]),digits = 3),
                                  ", ",round(mean(delta[,2]),digits = 3),")"))),
        side = 1, adj = 1, at = 3)
}