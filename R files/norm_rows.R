norm_rows <- function(network_list){
  normed_list = list()
  k = ncol(network_list[[1]]) 
  #temp_mat = matrix(nrow = k,ncol = k)
  for(i in 1:length(network_list)){
      normed_mat = normalize.rows(abs(network_list[[i]]), method = "manhattan")
      #temp_mat = abs(network_list[[i]])
      #normed_mat = normalize.rows(temp_mat, method = "maximum")
      #normed_mat = normalize.rows(temp_mat, method = "manhattan")
      rowSums(normed_mat)
      normed_list[[length(normed_list)+1]] = normed_mat
  }
  return(normed_list)
}
