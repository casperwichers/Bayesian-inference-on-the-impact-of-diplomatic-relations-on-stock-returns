plots_alphabeta <- function(alpha, beta){
  par(mfrow=c(4,3))
  for(i in 1:6){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha[,i], trim = 0.05)-0.1),
                      max(Trim(alpha[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_DGP[i], col="red")
  }
  for(i in 1:6){ # BETA PLOTS
    hist(Trim(beta[i,], trim = 0.02), freq = F, col = "blue",
         # breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.02) - 0.2),
         #              max(Trim(beta_results_fil[,i],trim = 0.02) + 0.2), by=0.1),
         breaks = 20,
         main = paste("Hist of MH for beta_tilde", i))
    abline(v = beta_tilde_DGP[i], col = "red")
  }
}