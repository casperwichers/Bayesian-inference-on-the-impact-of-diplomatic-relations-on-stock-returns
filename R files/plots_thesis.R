plots_country <- function(alpha, beta, delta, rho, sigma, country_list){
  par(mar=c(5,5,4,1)+.1)
  ####### ALPHA #######
  for(i in 1:length(country_list)){
    par(mfrow=c(6,3))
    ##### ALPHA #########
    plot(alpha[,i], type = "l", ylab = "Value", xlab = "",
         main = paste0("Plot of alpha results for ", country_list[i]), lwd=2,
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2)
    if(exists("alpha_DGP")){
      abline(h = alpha_DGP[i], col="red", lwd = 4)  
    }
    hist(alpha[,i], freq = F,
         breaks = seq(min(alpha[,i]-0.1),
                      max(alpha[,i]+0.1), by=0.05),
         main = paste0("Histogram of alpha results for ", country_list[i]),
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         xlab = "",
         col = "black")
    if(exists("alpha_DGP")){
      abline(v = alpha_DGP[i], col="red", lwd = 4)
    } else {
      abline(v = mean(alpha[,i]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(alpha[,i], main = paste0("ACF of alpha results for ", country_list[i]),
        ylab = "", cex.lab = 2.2, cex.axis = 2)
    ###### BETA #######
    
    ##### BETA_1
    plot(beta[,i], type = "l", ylab = "Value", cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         xlab = "",
         lwd=2, main = paste0("Plot of beta_oil results for ", country_list[i]))
    if(exists("beta_tilde_DGP")){
      abline(h = beta_tilde_DGP[i], col="red", lwd=4)  
    }
    hist(beta[,i], freq = F, cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
       # breaks = seq(min(beta[,i]-0.1),
       #              max(beta[,i]+0.1), by=0.05),
       breaks = 20,
       main = paste0("Histogram of beta_oil results for ", country_list[i]),
       xlab = "",
       col = "black")
    if(exists("beta_tilde_DGP")){
       abline(v = beta_tilde_DGP[i], col="red", lwd = 4)
    } else {
    abline(v = mean(beta[,i]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(beta[,i], cex.lab=2.2, cex.axis = 2,
        ylab = "",
        main = paste0("ACF of beta_oil results for ", country_list[i]))
    
    ##### BETA_2
    plot(beta[,(i + num_countries)], cex.main = 2.2, cex.lab = 2.2, cex.axis = 2,
         type = "l", ylab = "Value",lwd=2, xlab = "",
        main = paste0("Plot of beta_mat results for ", country_list[i]))
    if(exists("beta_tilde_DGP")){
      abline(h = beta_tilde_DGP[i + num_countries], col="red", lwd=4)  
    }
    hist(beta[,(i + num_countries)], freq = F,
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         # breaks = seq(min(beta[,(i + num_countries)]-0.1),
         #              max(beta[,(i + num_countries)]+0.1), by=0.05),
         breaks = 20,
         main = paste0("Histogram of beta_mat results for ", country_list[i]),
         xlab = "",
         col = "black")
    if(exists("beta_tilde_DGP")){
       abline(v = beta_tilde_DGP[i + num_countries], col="red", lwd = 4)
    } else {
    abline(v = mean(beta[,(i + num_countries)]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(beta[,(i + num_countries)], cex.lab = 2.2, cex.axis = 2,
        ylab = "",
        main = paste0("ACF of beta_mat results for ", country_list[i]))
    
    #### BETA_3
    plot(beta[,(i + 2*num_countries)], type = "l", ylab = "Value",lwd=2,
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         xlab = "",
         main = paste0("Plot of beta_GEPU results for ", country_list[i]))
    if(exists("beta_tilde_DGP")){
      abline(h = beta_tilde_DGP[i + 2*num_countries], col="red", lwd = 4)  
    }
    hist(beta[,(i + 2*num_countries)], freq = F,
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         # breaks = seq(min(beta[,(i + 2*num_countries)]-0.1),
         #              max(beta[,(i + 2*num_countries)]+0.1), by=0.05),
         breaks = 20,
         main = paste0("Histogram of beta_GEPU results for ", country_list[i]),
         xlab = "",
         col = "black")
    if(exists("beta_tilde_DGP")){
       abline(v = beta_tilde_DGP[i + 2*num_countries], col="red", lwd = 4) 
    } else {
      abline(v = mean(beta[,(i + 2*num_countries)]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(beta[,(i + 2*num_countries)], 
        cex.lab=2.2, cex.axis = 2,
        ylab = "",
        main = paste0("ACF of beta_GEPU results for ", country_list[i]))
  
    ##### RHO #######
    plot(rho[,i], type = "l", cex.main = 2.2, cex.lab= 2.2, cex.axis = 2, xlab = "",
         ylab = "Value",main = paste0("Plot of rho results for ", country_list[i]),lwd=2)
    if(exists("rho_DGP")){
      abline(h = rho_DGP[i], col="red", lwd=4)  
    }
    hist(rho[,i], col = "black", main = paste0("Histogram of rho results for ", country_list[i]),
         freq = F, 
         # breaks = seq(min(rho[,i]-0.1),
         #                        max(rho[,i]+0.1), by=0.01),
         breaks = 20,
         xlim = c(0,1),
         xlab = "",cex.main = 2.2, cex.lab= 2.2, cex.axis = 2)
    if(exists("rho_DGP")){
       abline(v = rho_DGP[i], col="red", lwd=4)
    } else {
    abline(v = mean(rho[,i]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(rho[,i], cex.lab = 2.2, cex.axis = 2, ylab = "",
        main = paste0("ACF of rho results for ", country_list[i]))
    
    ##### SIGMA #########
    plot(Sigma_results_fil[,i], type = "l", ylab = "Value", xlab = "",
         main = paste0("Plot of Sigma results for ", country_list[i]), lwd=2,
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2)
    if(exists("Sigma_DGP")){
      abline(h = Sigma_DGP[i,i], col="red", lwd = 4)  
    }
    hist(Sigma_results_fil[,i], freq = F,
         #breaks = seq(min(Sigma_results_fil[,i]-0.1),
         #              max(Sigma_results_fil[,i]+0.1), by=0.1),
         breaks = 20,
         main = paste0("Histogram of Sigma results for ", country_list[i]),
         cex.main = 2.2, cex.lab= 2.2, cex.axis = 2,
         xlab = "",
         col = "black")
    if(exists("Sigma_DGP")){
       abline(v = Sigma_DGP[i,i], col="red", lwd = 4) 
    } else {
      abline(v = mean(Sigma_results_fil[,i]), col="blue", lwd = 4, type = "l", lty = 2)
    }
    par(cex.main = 2.2)
    acf(Sigma_results_fil[,i], main = paste0("ACF of Sigma results for ", country_list[i]),
        ylab = "", cex.lab = 2.2, cex.axis = 2)
    
  }
}
 
plots_country(alpha_results_fil, beta_results_fil,
              delta_results_fil,rho_results_fil, Sigma_results_fil,
              country_list2)

