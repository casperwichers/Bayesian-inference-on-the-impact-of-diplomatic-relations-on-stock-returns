rm(list=ls())
cat("\014") 
library("plyr")

#read in required packages
require(readxl)
require(tidyverse)

#set the working directory from which the files will be read from
setwd("/Users/casper/Desktop/Bestanden")
#dataWeek <- read.table("20200505-icews-events.tab", header=T,sep = "\t", 
#                    fill = T, encoding = "automatic", quote = "\"")
            
 
# I will look at the observations during 2018:

#dataYear <- read.table("events.2018.20200427084805.tab", header=T,sep = "\t", 
#                    fill = T, encoding = "automatic", quote = "\"")

#dataYear_Fil = subset(dataYear, Source.Country != Target.Country)


# For all years:

file_list <- list.files(path="/Users/casper/Desktop/Bestanden")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- rep(0,20)

# had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <- read.table(file_list[i], header=T,sep = "\t", 
                          fill = T, encoding = "automatic", quote = "\"") 
  dataset <- rbind(dataset, temp_data)
}

dataYear_Fil = subset(dataset, Source.Country != Target.Country)

#Find Descriptives of filtered data:
range(dataYear_Fil[,"Intensity"]) # RESULT: -10 to 10
mean(dataYear_Fil[,"Intensity"])  # RESULT: 1.088588

sort(table(dataYear_Fil[,"Source.Country"]), decreasing = T)
# ### MOST OCCURING SOURCE COUTNRIES:
# USA:                  56288
# Russian Federation:   28163
# China:                20106
# India:                16331
# North Korea:          15262
# South Korea:          11551

# NaN Values: 23879

sort(table(dataYear_Fil[,"Target.Country"]), decreasing = T)
# ### MOST OCCURING TARGET COUTNRIES:
# USA:                  46699
# Russian Federation:   29789
# China:                19352
# North Korea:          18671
# Syria:                19352
# South Korea:          10680
# India:                10640

# NaN Values: 40026


# Since USA appears most, zoom in:
USA = subset(dataYear_Fil, dataYear_Fil[,"Source.Country"] == "United States")
sort(table(USA[,"Target.Country"]), decreasing = T)
# ### MOST OCCURING TARGET COUTNRIES WHERE SOURCE COUNTRY = USA:
# North Korea:          8647
# Russian Federation:   6976
# China:                3882
# South Korea:          2352
# Syria:                2204

# NaN Values: 3624

# Possibility: USA, Russia, Syria, North Korea, China

Russia = subset(dataYear_Fil, dataYear_Fil[,"Source.Country"] == "Russian Federation")
sort(table(Russia[,"Target.Country"]), decreasing = T)
# ### MOST OCCURING TARGET COUTNRIES WHERE SOURCE COUNTRY = RUSSIA:
# USA:                  5194
# Syria:                1915
# Turkey:               1795
# Ukraine:              1302
# China:                1049          
# North Korea:          959 

# NaN Values: 3624


China = subset(dataYear_Fil, dataYear_Fil[,"Source.Country"] == "China")
sort(table(China[,"Target.Country"]), decreasing = T)
# ### MOST OCCURING TARGET COUTNRIES WHERE SOURCE COUNTRY = China:
# USA:                  3733
# North Korea:          2052 
# Japan:                1252
# Russian Federation:   1028
# India:                881

# Syria:                33      <- This could be problematic

China_Syria = subset(China, China["Target.Country"] == "Syria")
mean(unlist(China_Syria["Intensity"]))
# However, the mean intensity of these
# observations is 2.79 which is much higher than average 

# NaN Values: 1457

NK = subset(dataYear_Fil, dataYear_Fil[,"Source.Country"] == "North Korea")
sort(table(NK[,"Target.Country"]), decreasing = T)
# ### MOST OCCURING TARGET COUTNRIES WHERE SOURCE COUNTRY = China:
# USA:                  6774
# South Korea:          4370 
# China:                1763
# Russian Federation:   688
         
# Syria:                89      <- This could be problematic

# NaN Values: 185

### CONCLUSION:
### IF POSSIBLE, MAYBE OPTIMAL TO USE USA, NORTH KOREA, CHINA AND RUSSIA;
### ALL THESE COUTNRIES HAVE A GOOD AMOUNT OF LINKS TO EACH OTHER

