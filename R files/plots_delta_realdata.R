plots_delta <- function(delta){
  par(mfrow=c(2,3))
  for(i in 1:2){
    plot(delta[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("δ", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(delta[,i])-0.1,
                max(delta[,i])+0.1))
    hist(delta[,i], freq = F,
         breaks = seq(min(delta[,i]) - 0.1,
                      max(delta[,i]) + 0.1, by=0.01),
         main = paste("Hist of MH values for",gsub(" ", "", paste("δ", i))),
         xlab = "Value",
         col = "blue")
  }
  # mtext(paste("n_iter =",niter, ", burn_in =", burn_in,
  #             ", start value = (0.25,0.25,0.25,0.25)",
  #             ", alpha =",gsub(" ", "", paste("(",nu[1], ",",nu[2],",", nu[3],",",nu[4],")"))),
  #       side = 3, line = -18, outer = TRUE)
}