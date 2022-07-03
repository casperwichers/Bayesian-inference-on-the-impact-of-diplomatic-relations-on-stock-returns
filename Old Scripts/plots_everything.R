plots_everything <- function(){
  par(mfrow=c(6,4))
  for(i in 1:4){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α^*", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_star_DGP[i], col="red")
  }
  for(i in 1:4){ # BETA PLOTS
    hist(Trim(beta_results_fil[,i], trim = 0.01), freq = F, col = "blue",
         breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.01) - 3),
                      max(Trim(beta_results_fil[,i],trim = 0.01) + 3), by=0.1),
         main = paste("Hist of sampled values of vec_beta^*", i))
    abline(v = vec_beta_star_DGP[i], col = "red")
  }
  for(i in 1:4){ # DELTA TRACE PLOTS
    plot(delta_results_fil[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results_fil[,i])-0.1,
                max(delta_results_fil[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
  }
  for(i in 1:4){ # DELTA HIST PLOTS
    hist(delta_results_fil[,i], freq = F,
         breaks = seq(min(delta_results_fil[,i]) - 0.1,
                      max(delta_results_fil[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i)), "nu = ", nu[1]),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
    #acf(delta_results_fil[,i], lag.max = 10, type = "correlation", plot = TRUE)
  }
  for(i in 1:6){ # RHO PLOTS
    hist(rho_results_fil[,i], freq = F, col = "blue",
         breaks = seq(min(rho_results_fil[,i]) - 0.1, 
                      max(rho_results_fil[,i]) + 0.1, by=0.005),
         main = paste0("Histogram of results of rho", i))
    abline(v = rho_DGP[i], col="red")
  }
  mtext(paste("# Iterations =", nreps, ", # Burned-in samples =", burn_in,
              "\n start value = (0.25,0.25,0.25,0.25)",
              ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")")),
              "\n delta means = ", 
              gsub(" ", "", paste("(",round(mean(delta_results_fil[,1]),digits = 3),
                             ", ",round(mean(delta_results_fil[,2]),digits = 3),
                             ", ",round(mean(delta_results_fil[,3]),digits = 3),
                             ", ",round(mean(delta_results_fil[,4]),digits = 3),")"))),
        at = 1,1)
}