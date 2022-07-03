plots_alpha <- function(alpha){
  par(mfrow=c(6,2))
  for(i in 1:ncol(alpha)){
    plot(alpha[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("α^*", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(alpha[,i])-0.5,
                max(alpha[,i])+0.5))
    abline(h = alpha_star_DGP[i], col="red")
    hist(Trim(alpha[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(alpha[,i], trim = 0.05)-0.1),
                      max(Trim(alpha[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("α", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = alpha_star_DGP[i], col="red")
  }
}