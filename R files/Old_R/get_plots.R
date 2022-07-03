plot_counts <- function(data, dates){
  counts = vector()
  for(i in 1:(length(dates)-1)){
    selection = data[(data$Event.Date >= dates[i]) & 
                       (data$Event.Date < dates[i+1]),]
    counts[length(counts)+1] = nrow(selection)
  }
  plot(x = head(dates,-1),
       y = counts,
       xlab = "Date",
       ylab = "Count",
       main = "Monthly frequency of events")
}

plot_counts(data, dates)

