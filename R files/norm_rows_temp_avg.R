norm_rows_temp_avg <- function(network_list){
  normed_list = list()
  k = ncol(network_list[[1]]) 
  temp_mat = matrix(0, nrow = k,ncol = k)
  for(t in 1:length(network_list)){
    temp_mat = temp_mat + network_list[[t]]
    }
  normed_mat = normalize.rows(abs(temp_mat), method = "manhattan")
  for(t in 1:length(network_list)){
    normed_list[[t]] = normed_mat
  }
  return(normed_list)
}
