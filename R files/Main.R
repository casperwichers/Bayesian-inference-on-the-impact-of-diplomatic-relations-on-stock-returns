rm(list=ls())
cat("\014") 
library("plyr")
library(beepr)
library(Matrix)
library(wordspace)

# Read in Files
data = read.csv('/Users/casper/Desktop/Scriptie/csv files/data_clean.csv')
data = subset(data, Quad.Class.Hoff != "NaN")

#start_date = min(data$Event.Date)
#end_date = max(data$Event.Date)

months.list <- c("01","02","03","04",
                 "05","06","07","08",
                 "09","10","11","12")
years.list = c(2000:2016)
df<-expand.grid(year=years.list,month=months.list)
df$Date=as.Date(paste0(df$year,"-",df$month,"-01"))
df<-df[order(df$Date),]
dates = df$Date
dates = head(dates,-11)
T = length(dates)-1
start_date = min(dates)
end_date = max(dates)

matrices = c("m_min_mat", "m_plus_mat", "v_min_mat", "v_plus_mat")

###### Choose countries ######### 
country_list = list("United Kingdom",
                    "China",
                    "India",
                    "Russian Federation",
                    "Japan",
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
                    "United States")

selected = subset(data, Source.Country %in% country_list)
data = subset(selected, Target.Country %in% country_list)

write.csv(x = data, 
          "/Users/casper/Desktop/Scriptie/csv files/data_countryselection.csv",
          row.names = FALSE)

#############################
num_countries = length(country_list)

m_min_list = list()
m_plus_list = list()
v_min_list = list()
v_plus_list = list()

counts = vector()

for(i in 1:(length(dates)-1)){
  for(j in matrices){
    assign(j, matrix(data = rep(0),
                     nrow=length(country_list),
                     ncol=length(country_list),
                     byrow=FALSE,
                     dimnames = list(c(unlist(country_list)),
                                     (c(unlist(country_list))))))
  }
  selection = data[(data$Event.Date >= dates[i]) & 
                     (data$Event.Date < dates[i+1]),]
  counts[length(counts)+1] = nrow(selection)
  
  if(nrow(selection > 0)){
    for(k in 1:nrow(selection)){
       source = selection[k,"Source.Country"]
       target = selection[k, "Target.Country"]
       class = selection[k, "Quad.Class.Hoff"]
       intensity = selection[k,"Intensity"]

       if(class == "v+"){
         v_plus_mat[source,target] = v_plus_mat[source,target] + intensity
       } else if(class == "v-"){
         v_min_mat[source,target] = v_min_mat[source,target] + intensity
       } else if(class == "m+"){
         m_plus_mat[source,target] = m_plus_mat[source,target] + intensity
       } else if(class == "m-"){
         m_min_mat[source,target] = m_min_mat[source,target] + intensity
       }
    }
  }
  m_min_list[[length(m_min_list)+1]] = m_min_mat
  m_plus_list[[length(m_plus_list)+1]] = m_plus_mat
  v_min_list[[length(v_min_list)+1]] = v_min_mat
  v_plus_list[[length(v_plus_list)+1]] = v_plus_mat
}

##### MAX NORMS ######
# max_norm_rows <- function(network){
#   T = length(network)
#   num_countries = nrow(network[[1]])
#   colsums = rep(0,num_countries)
#   for(t in 1:T){
#     for(i in 1:num_countries){
#       print(unname(colSums(network[[t]])[i]))
#       if(abs(unname(colSums(network[[t]])[i])) > abs(colsums[i])){
#         colsums[i] = colSums(network[[t]])[i]
#         #print("yoyoyo")
#       }
#     }
#   }
#   for(t in 1:T){
#     for(i in 1:num_countries){
#       network[[t]][,i] = network[[t]][,i] / colsums[i] 
#     }
#   }
#   print(colsums)
#   return(network)
# }

# m_min_normed = max_norm_rows(m_min_list)
# m_plus_normed = max_norm_rows(m_plus_list)
# v_min_normed = max_norm_rows(v_min_list)
# v_plus_normed = max_norm_rows(v_plus_list)
###################

beep(sound = 1)
#rm(selected)
rm(df,selection,class,i,intensity,months.list,source,target,years.list,j,k)
rm(selected,matrices)
rm(m_min_mat,m_plus_mat,v_min_mat,v_plus_mat)


