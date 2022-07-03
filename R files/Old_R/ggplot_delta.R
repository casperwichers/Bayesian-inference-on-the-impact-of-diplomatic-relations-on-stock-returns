library(ggplot2)
# Basic line plot with points
ggplot_delta <- function(delta){
  df = as.data.frame(delta)
  #df$index = rownames(df)
  head(df)
  
  #qplot(rownames(df),df$V1, geom="auto")
  ggplot(data=df, aes(x=rownames(df), y=V1, group=1)) +
    geom_line()+
    geom_point()

}


ggplot_delta(delta_results)