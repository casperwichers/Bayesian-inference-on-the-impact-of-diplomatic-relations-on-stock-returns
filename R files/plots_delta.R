plots_delta <- function(delta){
  par(mfrow=c(ncol(delta),3))
  for(i in 1:ncol(delta)){
    plot(delta[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "black",
         cex.main = 1.5,
         cex.lab = 1.5,
         ylim=c(min(delta[,i])-0.1,
                max(delta[,i])+0.1))
    if(exists("Delta_DGP")){
      abline(h = delta_DGP[i], col="red")
    }
    hist(delta[,i], freq = F,
         breaks = seq(min(delta[,i]) - 0.1,
                      max(delta[,i]) + 0.1, by=0.02),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
         xlab = "Value",
         cex.main = 1.5,
         cex.lab = 1.5,
         col = "black")
    if(exists("Delta_DGP")){
      abline(v = delta_DGP[i], col="red")  
    }
    temp = acf(delta[,i], lag.max = 15, type = "correlation", plot=FALSE)
    plot(temp, 
        cex.main = 2,
        cex.lab = 3,
        main = paste("ACF for",gsub(" ", "", paste("δ", i))))
  }
}
