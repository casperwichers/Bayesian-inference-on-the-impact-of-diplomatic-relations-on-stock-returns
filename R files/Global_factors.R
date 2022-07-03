#rm(list=ls())
#cat("\014") 
library("plyr")
library(dplyr)
library(xts)
require(data.table)
require(zoo)

#start_date = min(dates)
#end_date = max(dates)

globalfactors <- function(ticker, start_date,end_date){
  link = gsub(" ", "",
              paste('/Users/casper/Desktop/Scriptie/Global_factors_data/',
                    ticker,'.csv'))
  if(ticker == 'POILBREUSDM'){
    ticker = 'Global Brent Crude Oil price'
  } else if(ticker == 'PINDUINDEXM'){
    ticker = 'Global Industrial Materials Index'
  } else if(ticker == 'GEPUPPP'){
    ticker = 'PPP-Adjusted GDP'
  }
  data = read.table(link, sep = ",", header = TRUE)
  data[,1] = as.Date(data[,1])
  temp = diff(log(data[,2]))
  data = head(data, -1)
  data[,2] = temp
  data = subset(data, DATE >= start_date & DATE < end_date)
  #print(data)
  # plot(x = data[,1], y = data[,2],
  #      xlab = 'Date', ylab = 'Price (in USD)',
  #      type="l", lwd = "1",
  #      main = paste0("Plot of ",ticker," data"))
  return(data)
}

#par(mfrow=c(1,3), cex=0.7, mai=c(0.65,0.7,0.3,0.3))
oil = globalfactors('POILBREUSDM',start_date, end_date)
mat_idx = globalfactors('PINDUINDEXM',start_date, end_date)
GDP = globalfactors('GEPUPPP', start_date, end_date)

global_factors = merge(oil, c(mat_idx, GDP), by = "DATE", all.x = TRUE)
global_factors <- global_factors[-4]
names(global_factors) <- c("Date","Oil","Materials index","GEPU")
facts <- xts(global_factors[,-1], order.by=as.Date(global_factors[,1], "%Y/%m/%d"))
rownames(global_factors) <- global_factors$Date; global_factors<-global_factors[,-1]

if(sum(is.nan(as.matrix(global_factors))) == 0){
  write.csv(x = global_factors, 
            "/Users/casper/Desktop/Scriptie/csv files/global_factors.csv",
            row.names = TRUE)
}

plot.xts(facts[,1:ncol(facts)],
     main = "Monthly log-differences of the global factors",
     lwd = 3,
     cex = 1.2)
addLegend("topright", on=1, 
          lty=c(1, 1), lwd=c(2, 1),
          box.lwd = par("lwd"), 
          box.lty = par("lty"), box.col = par("fg"),
          cex = 1
)

rm(GDP, mat_idx, oil)
