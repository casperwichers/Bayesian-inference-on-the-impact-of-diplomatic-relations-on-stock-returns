plots_alphabetadelta <- function(){
  par(mfrow=c(6,3))
  for(i in 1:5){ # ALPHA PLOTS
    # plot(alpha_results_fil[,i], type = "l", lwd = 3,
    #      main = paste("Plot of MH values for",gsub(" ", "", paste("α", i))),
    #      ylab = "Value",
    #      col = "blue",
    #      ylim=c(min(alpha_results[,i])-0.2,
    #             max(alpha_results[,i])+0.2))
    # abline(h = alpha_DGP[i], col="red")
    hist(Trim(alpha_results_fil[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha_results_fil[,i], trim = 0.05)-0.1),
                      max(Trim(alpha_results_fil[,i], trim = 0.05)+0.1), by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α^*_", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_star_DGP[i], col="red")
  }
  for(i in 1:5){ # BETA PLOTS
    hist(Trim(beta_results_fil[,i], trim = 0.05), freq = F, col = "blue",
         breaks = seq(min(Trim(beta_results_fil[,i], trim = 0.05) - 0.2),
                      max(Trim(beta_results_fil[,i],trim = 0.05) + 0.2), by=0.1),
         main = paste("Hist of sampled values of vec(beta_star)", i))
    abline(v = vec_beta_star_DGP[i], col = "red")
  }
  for(i in 1:4){
    plot(delta_results[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta_results[,i])-0.1,
                max(delta_results[,i])+0.1))
    abline(h = delta_DGP[i], col="red")
    hist(delta_results_fil[,i], freq = F,
         breaks = seq(min(delta_results_fil[,i]-0.05),
                      max(delta_results_fil[,i])+0.05, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
         xlab = "Value",
         col = "blue")
    abline(v = delta_DGP[i], col="red")
  }
}