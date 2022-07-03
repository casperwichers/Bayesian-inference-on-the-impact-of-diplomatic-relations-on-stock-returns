plot_counts <- function(counts, dates){
  # for(i in 1:(length(dates)-1)){
  #   selection = data[(data$Event.Date >= dates[i]) & 
  #                      (data$Event.Date < dates[i+1]),]
  #   counts[length(counts)+1] = nrow(selection)
  plot(x = head(dates,-1),
       y = counts,
       xlab = "Date",
       lwd = "2",
       pch = 19,
       ylab = "Monthly count",
       cex.lab=2,
       cex.axis=1.5,
       cex.main = 2,
       main = "Monthly frequency of events")
}

plot_counts(counts, dates)



library(extrafont)
font_install('fontcm')
loadfonts()

check = head(dates,-1)
df = as.data.frame(cbind(check, counts), row.names = check)
#df$check1 = as.Date(origin = df$check)
head(df)
df$idu = row.names(df)
library(ggplot2)
a = ggplot(data = df,
       aes(x=check, y=counts)) + geom_point(size=3) +
  #ggtitle("Monthly count of observations during the period 2000-01-01/2015/12/31") +
  xlab("Date") + ylab("Frequency") +
  theme_bw() + 
  theme(axis.text=element_text(size=16, family= "CM Roman"), axis.title = element_text(size=20)) + 
  scale_x_date(date_breaks = "36 months", date_labels =  "%Y-%m-%d") 
print(a)

