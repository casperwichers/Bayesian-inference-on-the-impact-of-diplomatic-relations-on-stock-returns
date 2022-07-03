get_dates <- function(start_year,end_year){
  months.list <- c("01","02","03","04",
                   "05","06","07","08",
                   "09","10","11","12")
  years.list <- c(start_year:(end_year+1))
  df<-expand.grid(year=years.list,month=months.list)
  df$Date=as.Date(paste0(df$year,"-",df$month,"-01"))
  df<-df[order(df$Date),]
  dates = df$Date
  dates = head(dates,-11)
  start_date = min(dates)
  end_date = max(dates)
  T = length(dates)-1
  return(list(dates, start_date, end_date, T))
}

