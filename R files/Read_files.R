rm(list=ls())
cat("\014") 
library("plyr")
library(beepr)

#read in required packages
require(readxl)
require(tidyverse)

library(lubridate)
library(dplyr)

################################################################################
# READ IN TAB FILES:
#set the working directory from which the files will be read from
setwd("/Users/casper/Desktop/TVP-SAR/Bestanden")
#create a list of the files from your target directory
file_list <- list.files(path="/Users/casper/Desktop/TVP-SAR/Bestanden")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- rep(0,20)

# had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  print(i)
  temp_data <- read.table(file_list[i], 
                          header=T,
                          sep = "\t", 
                          fill = T, 
                          encoding = "automatic", 
                          quote = "")
  
  # Remove quotes if any left
  character_cols = which(sapply(temp_data, class) == 'character')
  for(i in 1:length(character_cols)) {
    a = character_cols[i]
    temp_data[,a] = gsub("\"", "", temp_data[,a])
  }   
  #Remove zeros before CAMEO:
  temp_data[,"CAMEO.Code"] <- gsub('^0','',temp_data[,"CAMEO.Code"])
  dataset <- rbind(dataset, temp_data)
}

#read in types 
types = read.csv('/Users/casper/Desktop/Scriptie/csv files/events_count_cameo.csv',
                 header = TRUE,sep=";")
types = types[,4:5]

################################################################################

# Choose countries
#country_list = list(#"United Kingdom",
                    #"China",
                    #"India",
                    #"Russian Federation",
                    #"Japan",
                    #"Iran",
                    #"Australia",
                    #"Israel",
                    #"France",
                    #"Germany",
                    #"Brazil",
                    #"Canada",
                    #"Turkey",
                    #"Taiwan",
                    #"South Korea",
                    #"United States")
# 
# selected = subset(dataset, Source.Country %in% country_list)
# dataset = subset(selected, Target.Country %in% country_list)

################################################################################

# Remove columns and non-international observations
#data_filtered = subset(dataset, Source.Country != Target.Country)
data_filtered = dataset
data_filtered = subset(data_filtered, select = -c(Event.ID, Story.ID, Latitude,Longitude,
                                                  City,District,Province,
                                                  Sentence.Number,Publisher,City,District,
                                                  Province,Source.Sectors,Country,
                                                  Target.Sectors))

################################################################################

# Order columns
#data = data_filtered[1:1000,]
joined = merge(data_filtered,types,by.y='CAMEO.Code', all.x = T)
data = joined[c(2,1,5,9,6,3,4,7,8)]
# Sort on date
data$Event.Date <- as.Date(data$Event.Date, format = "%Y-%m-%d")
#class(data$Event.Date)
data <- data[order(data$Event.Date),] 

################################################################################

sort(table(data[,"Source.Country"]), decreasing = T)
sort(table(data[,"Target.Country"]), decreasing = T)

################################################################################

# Export to CSV file
write.csv(x = data, 
          "/Users/casper/Desktop/Scriptie/csv files/data_clean.csv",
          row.names = FALSE)

rm(data_filtered,dataset,joined,temp_data,types,a,character_cols,file_list,i)
beep(sound = 1)
