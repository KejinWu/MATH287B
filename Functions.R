# Define the similarity function

similarity = function(x1 = NULL, x2 = NULL, sigma=1) {
  similarity_value = exp(- norm(as.matrix(x1-x2), type="F")^2/(2*sigma^2))
}

# Build Similarity graph
Build_similarity_graph = function(data = NULL) {
  if(is.matrix(data)){
  N = nrow(data)
  S = matrix(nrow = N, ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] = similarity(data[i,], data[j,])
    }
  }
return("Similarity_graph" = S)
  }
  
  if(is.vector(data)){
    N = length(data)
    S = matrix(nrow = N, ncol=N)
    for(i in 1:N) {
      for(j in 1:N) {
        S[i,j] = similarity(x1 = data[i], x2 = data[j])
      }
    }
    return("Similarity_graph" = S)
    
  }
}

# Build W matrix

Build_W = function(Similarity_graph = NULL, k_near=NULL, mutual_k_near = NULL,epsilon = NULL, Full_graph = FALSE) {
  N = ncol(Similarity_graph)
  
  # Build fully connected graph
  if (Full_graph == TRUE) {  
    W = Similarity_graph
  } 

  # Build k-near graph
if (is.null(k_near) == FALSE) {
    W = matrix(0,nrow = N, ncol=N)
    for(i in 1:N) { 
      Order_similarity = sort(Similarity_graph[i,], decreasing=TRUE)[1:k_near]
      for (k_near_value in Order_similarity) {
        j = which(Similarity_graph[i,] == k_near_value)
        W[i,j] = Similarity_graph[i,j]
        W[j,i] = W[i,j]
      }
    }
}
  # Build mutual k-near graph
  if (is.null(mutual_k_near) == FALSE) {
    W = matrix(0,nrow = N, ncol=N)
    connected_info = vector("list",length = N)
    for(i in 1:N) { 
      Order_similarity = sort(Similarity_graph[i,], decreasing=TRUE)[1:mutual_k_near]
      connected_points = c()
      for (k_near_value in Order_similarity) {
        j = which(Similarity_graph[i,] == k_near_value)
        connected_points = c(connected_points,j)
      }
      connected_info[[i]] = connected_points
    }
    
    for (i in 1:N) {
      for (j in connected_info[[i]]) {
        if(i%in%connected_info[[j]]){
          W[i,j] = Similarity_graph[i,j]
        }
      }
    }
  }
  
  # Build k-near graph
  if (is.null(epsilon) == FALSE) {
    W = matrix(0,nrow = N, ncol=N)
    for(i in 1:N) { 
      for (j in (i):N) {
        if(Similarity_graph[i,j]>(1-epsilon)){
        W[i,j] = Similarity_graph[i,j]
        W[j,i] = W[i,j]
        }
      }
    }
  }
  
  return("W" = W)
}


# Perform spectral clustering (Unnormalized & Normalized )

Spectral_clustering = function(data = NULL,sigma = 1,k_near=NULL, mutual_k_near = NULL,epsilon = NULL, Full_graph = FALSE,Graph_laplacian = "Normalized_rw",Number_of_cluster = NULL){
  S = Build_similarity_graph(data)
  W = Build_W(Similarity_graph = S,epsilon =  epsilon,k_near = k_near, mutual_k_near = mutual_k_near,Full_graph = Full_graph)
  N = nrow(W)
  D = diag(apply(W, 1, sum))
  if (Graph_laplacian == "Unnormalized"){
  L = D - W
  eigen_vectors = eigen(L,symmetric=TRUE)$vectors
  eigen_values = eigen(L,symmetric=TRUE)$values
  selected_vectors = eigen_vectors[,(ncol(eigen_vectors)-Number_of_cluster+1):ncol(eigen_vectors)]
  km_cluser = kmeans(selected_vectors, centers=Number_of_cluster, nstart=1)
  plot(data, col=km_cluser$cluster,ylab = "Y value",xlab = "X value",main = "Clustering results \n unnormalized graph laplacian")
  }

  if (Graph_laplacian == "Normalized_rw"){
    L = diag(1,N) - mpower(D,-1)%*%W
    eigen_vectors = eigen(L,symmetric=TRUE)$vectors
    eigen_values = eigen(L,symmetric=TRUE)$values
    selected_vectors = eigen_vectors[,(ncol(eigen_vectors)-Number_of_cluster+1):ncol(eigen_vectors)]
    km_cluser = kmeans(selected_vectors, centers=Number_of_cluster, nstart=1)
    plot(data, col=km_cluser$cluster,ylab = "Y value",xlab = "X value",main = "Clustering results \n normalized graph laplacian")
  }
  return(list("eigen_vectors" = eigen_vectors,"eigen_values" = eigen_values,"selected_vectors" = selected_vectors,"D" = D))

}


# Unnormalized spectral clustering method with more information

Spectral_clustering_moreeigenvectors = function(data = NULL,more = 2,sigma = 1,k_near=NULL, mutual_k_near = NULL,epsilon = NULL, Full_graph = FALSE,Graph_laplacian = "Normalized_rw",Number_of_cluster = NULL){
  S = Build_similarity_graph(data)
  W = Build_W(Similarity_graph = S,epsilon =  epsilon,k_near = k_near, mutual_k_near = mutual_k_near,Full_graph = Full_graph)
  N = nrow(W)
  D = diag(apply(W, 1, sum))
  if (Graph_laplacian == "Unnormalized"){
    L = D - W
    eigen_vectors = eigen(L,symmetric=TRUE)$vectors
    eigen_values = eigen(L,symmetric=TRUE)$values
    selected_vectors = eigen_vectors[,(ncol(eigen_vectors)-Number_of_cluster+1-more):ncol(eigen_vectors)]
    km_cluser = kmeans(selected_vectors, centers=Number_of_cluster, nstart=1)
    plot(data, col=km_cluser$cluster,ylab = "Y value",xlab = "X value",main = paste("Use",more+2,"eigenvectors"),cex.main = 0.8)
  }
  
  if (Graph_laplacian == "Normalized_rw"){
    L = diag(1,N) - mpower(D,-1)%*%W
    eigen_vectors = eigen(L,symmetric=TRUE)$vectors
    eigen_values = eigen(L,symmetric=TRUE)$values
    selected_vectors = eigen_vectors[,(ncol(eigen_vectors)-Number_of_cluster+1):ncol(eigen_vectors)]
    km_cluser = kmeans(selected_vectors, centers=Number_of_cluster, nstart=1)
    plot(data, col=km_cluser$cluster,ylab = "Y value",xlab = "X value",main = "Clustering results \n normalized graph laplacian")
  }
  return(list("eigen_vectors" = eigen_vectors,"eigen_values" = eigen_values,"selected_vectors" = selected_vectors,"D" = D))
  
}

















