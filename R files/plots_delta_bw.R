plots_delta_bw <- function(delta){
  par(mfrow=c(2,3))
  #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
  par(mar=c(5,5,4,1)+.1)
  
  for(i in 1:2){
    plot(delta[,i], type = "l", lwd = 2,
         main = paste("Plot of MH values for delta", i),
         xlab = "",
         ylab = "Value",
         col = "black",
         cex.main = 2.5,
         cex.lab = 2.5,
         cex.axis = 2.2,
         ylim=c(0,1))
    if(exists("delta_DGP")){
      abline(h = delta_DGP[i], col="red", lwd = 4)
    }
    hist(delta[,i], plot = TRUE,
         xlab = "",
         col = "black",
         cex.main = 2.5,
         cex.lab = 2.5,
         xlim = c(0,1),
         cex.axis = 2.2,
         freq = F,
         main = paste("Histogram of results for delta",i)
    )
    if(exists("delta_DGP")){
      abline(v = delta_DGP[i], col="red", lwd = 4)  
    }
    par(cex.main = 2.5)
    acf(delta[,i], lag.max = 15, type = "correlation",
        cex.lab = 2.5,
        cex.axis = 2.2,
        ylab = "",
        main = paste("ACF of values of delta",i),
        plot=TRUE)
  }
  
  ###################
  # plot(delta[,1], type = "l", lwd = 2,
  #      main = "Plot of MH values for delta 1",
  #      ylab = "Value",
  #      col = "black",
  #      cex.main = 2.5,
  #      cex.lab = 2.5,
  #      cex.axis = 2,
  #      ylim=c(min(delta[,1])-0.1,
  #             max(delta[,1])+0.1))
  # if(exists("delta_DGP")){
  #   abline(h = delta_DGP[1], col="red",lwd=3)
  # }
  # plot(delta[,2], type = "l", lwd = 2,
  #      main = "Plot of MH values for delta 2",
  #      ylab = "Value",
  #      col = "black",
  #      cex.main = 2.5,
  #      cex.lab = 2.5,
  #      cex.axis = 2,
  #      ylim=c(min(delta[,2])-0.1,
  #             max(delta[,2])+0.1))
  # if(exists("delta_DGP")){
  #   abline(h = delta_DGP[2], col="red",lwd=3)
  # }
  # 
  # hist(delta[,1], plot = TRUE,
  #              xlab = "Value",
  #               col = "black",
  #               cex.main = 2.5,
  #               cex.lab = 2.5,
  #               cex.axis = 2,
  #               freq = F,
  #              main = "Histogram of results for delta 1"
  #              )
  # if(exists("delta_DGP")){
  #   abline(v = delta_DGP[1], col="red", lwd = 3)  
  # }
  # hist(delta[,2],  plot = TRUE,
  #       xlab = "Value",
  #       col = "black",
  #       cex.main = 2.5,
  #       cex.lab = 2.5,
  #       cex.axis = 2,
  #       freq = F,
  #      main = "Histogram of results for delta 2")
  # if(exists("delta_DGP")){
  #   abline(v = delta_DGP[2], col="red", lwd = 3)  
  # }
  # # breaks = seq(min(delta[,i]) - 0.1,
  # #              max(delta[,i]) + 0.1, by=0.02),
  # # main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
  # # xlab = "Value",
  # # cex.main = 1.5,
  # # cex.lab = 1.5,
  # # col = "dodgerblue")
  # 
  # par(cex.main=2)
  # acf(delta[,1], lag.max = 15, type = "correlation",
  #     #cex.main = 25,
  #     cex.lab = 2,
  #     cex.axis = 2,
  #     cex = 3,
  #     main = "ACF of values of delta",
  #     plot=TRUE)
  ######################
}

plots_delta_bw(delta_results_fil)
