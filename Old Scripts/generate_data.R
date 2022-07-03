#library(rowNorms)
library(wordspace)

#############################################################################
n = length(m_min_list)

# GENERATE DATA
v_plus_rows = list()
v_min_rows = list()
m_plus_rows = list()
m_min_rows = list()

for(i in 1:length(country_list)){
  country_matrix = matrix(ncol = length(country_list))
  v_plus_rows[[length(v_plus_rows)+1]] = country_matrix
  v_min_rows[[length(v_min_rows)+1]] = country_matrix
  m_plus_rows[[length(m_plus_rows)+1]] = country_matrix
  m_min_rows[[length(m_min_rows)+1]] = country_matrix
}

for(i in 1:n){
  for(j in 1:length(country_list)){
    v_plus_rows[[j]] = rbind(v_plus_rows[[j]], v_plus_list[[i]][j,])
    v_min_rows[[j]] = rbind(v_min_rows[[j]], v_min_list[[i]][j,])
    m_plus_rows[[j]] = rbind(m_plus_rows[[j]], m_plus_list[[i]][j,])
    m_min_rows[[j]] = rbind(m_min_rows[[j]], m_min_list[[i]][j,])
  }
}
# Remove first rows that contain only Nan
v_plus_rows = lapply(v_plus_rows, function(X) X[-1,])
v_min_rows = lapply(v_min_rows, function(X) X[-1,])
m_plus_rows = lapply(m_plus_rows, function(X) X[-1,])
m_min_rows = lapply(m_min_rows, function(X) X[-1,])

# Remove rows with a NaN
for(j in 1:length(country_list)){
    v_plus_rows[[j]] = na.omit(v_plus_rows[[j]])
    v_min_rows[[j]] = na.omit(v_min_rows[[j]])
    m_plus_rows[[j]] = na.omit(m_plus_rows[[j]])
    m_min_rows[[j]] = na.omit(m_min_rows[[j]])
}

sim_data <- function(rows_list, num_sims){
  #n = length
  sims = list()
  k = ncol(rows_list[[1]]) 
  #n = nrow(rows_list[[1]])
  temp_mat = matrix(nrow = k,ncol = k)
  for(i in 1:num_sims){
    for(j in 1:k){
      temp_mat[j,] = rows_list[[j]][max(1,round(runif(1)*nrow(rows_list[[j]]))),]
    }
    sims[[length(sims)+1]] = temp_mat
  }
  return(sims)
}

rm(country_matrix,i,j,n)

# v_plus_sims = sim_data(v_plus_rows,1000)
# v_min_sims = sim_data(v_min_rows,1000)
# m_plus_sims = sim_data(m_plus_rows,1000)
# m_min_sims = sim_data(m_min_rows,1000)

