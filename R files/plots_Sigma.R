plots_Sigma <- function(Sigma){
  par(mfrow=c(6,2))
  for(i in 1:ncol(Sigma)){
    print(Sigma_DGP[i,i])
    plot(Sigma[,i], type = "l", lwd = 3,
         main = paste("Plot of MH values for",gsub(" ", "", paste("σ^2_", i))),
         ylab = "Value",
         col = "blue",
         ylim=c(min(Sigma[,i])-0.5,
                max(Sigma[,i])+0.5))
    abline(h = Sigma_DGP[i,i], col="red")
    hist(Trim(Sigma[,i], trim = 0.05), freq = F,
         breaks = seq(min(Trim(Sigma[,i], trim = 0.05)-0.1),
                      max(Trim(Sigma[,i], trim = 0.05)+0.1), by=0.05),
         main = paste("Hist of MH values for",gsub(" ", "", paste("σ^2_", i)),"trimmed by 0.01"),
         xlab = "Value",
         col = "blue")
    abline(v = Sigma_DGP[i,i], col="red")
  }
}