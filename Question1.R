# install.packages('igraph')
# install.packages('cccd')
main <- function() {
  library('cccd')
  library('igraph')
  
  datafile <- read.csv('../Datasets/xAPI-edu-Data.csv')
  pearsons <- question_1a(datafile)
  question_1b(datafile)
  datafile <- read.csv("../Datasets/AlzheimersDisease.csv", header = TRUE)
  # generate_matrixes(datafile)
  get_relative_neighbourhoods()
}
 
get_relative_neighbourhoods <- function(){
  samples_matrix <- read.csv("../Datasets/samples_distance_matrix.csv", header = TRUE)
  proteins_matrix <- read.csv("../Datasets/proteins_distance_matrix.csv")
  rownames(samples_matrix) <- samples_matrix[,1]
  rownames(proteins_matrix) <- proteins_matrix[,1]
  proteins_matrix[,1] <- NULL
  samples_matrix[,1] <- NULL
  relative_neighbourhoods(samples_matrix,"samples")
  relative_neighbourhoods(proteins_matrix,"proteins")
}
relative_neighbourhoods <- function(matrix,title)
{
  print("hello")
  rng_graph <- rng(dx = matrix, algorithm = 'kd_tree', open = TRUE)
  rng_graph <- as.undirected(rng_graph)
  #adds labels
  V(rng_graph)$label <- colnames(matrix)
  
  E(rng_graph)$weight <- apply(get.edges(rng_graph, 1:gsize(rng_graph)),1,function(x)matrix[x[1],x[2]])
  E(rng_graph)$label <- E(rng_graph)$weight
  name <- paste("Results/",title,"RNG.gml", sep ="")
  write_graph(rng_graph, name, format = 'gml')
}
generate_matrixes <- function(datafile)
{
  rownames(datafile) <- datafile[,1]
  datafile[,1] <- NULL
  samples_matrix <- get_distance_samples(datafile)
  proteins_matrix <- get_distance_proteins(datafile)
  write.csv(samples_matrix,file ="../Datasets/samples_distance_matrix.csv")
  write.csv(proteins_matrix,file ="../Datasets/proteins_distance_matrix.csv")
  
  samples_graph <- mst(graph_from_adjacency_matrix(samples_matrix, mode = "undirected",weighted = TRUE))
  proteins_graph <- mst(graph_from_adjacency_matrix(proteins_matrix, mode = "undirected",weighted = TRUE))
  
  V(samples_graph)$label <- colnames(samples_matrix)
  E(samples_graph)$label <- round(E(samples_graph)$weight, 3)
  V(proteins_graph)$label <- colnames(proteins_matrix)
  E(proteins_graph)$label <- round(E(proteins_graph)$weight, 3)


  write_graph(samples_graph, 'Results/samplesMST.gml',format = 'gml')
  write_graph(proteins_graph, 'Results/proteinsMST.gml',format = 'gml')
  print("hello")
}
get_distance_proteins <- function(datafile)
{
  proteins <- matrix(nrow = nrow(datafile), ncol = nrow(datafile))
  colnames(proteins) <- rownames(datafile)
  rownames(proteins) <- rownames(datafile)
  i <- 1
  a <- 1
  print(proteins)
  while(i < nrow(proteins)+1){
    while(a < ncol(proteins)+1)
    {
      proteins[i,a] <- euclidean_dist(datafile[i,],datafile[a,])
      a = a + 1
    }
    i = i + 1
    a = 1
  }
  return(proteins)
}
get_distance_samples <- function(datafile)
{
  
  samples <- matrix(nrow = ncol(datafile), ncol = ncol(datafile))
  colnames(samples) <- colnames(datafile)
  rownames(samples) <- colnames(datafile)
  i <- 1
  a <- 1
  while(i < nrow(samples)+1){
    while(a < ncol(samples)+1)
    {
      samples[i,a] <- euclidean_dist(datafile[,i],datafile[,a])
      a = a + 1
    }
    i = i + 1
    a = 1
  }
  return(samples)
}



question_1a <- function(datafile){
  #function takes a dataframe as input and uses the cor method which calculates the pearsons correlation for all the variable pairs.
  numeric_columns <- datafile[ names(datafile) %in% c("raisedhands","VisITedResources","AnnouncementsView","Discussion")]
  pearsons <- cor(numeric_columns, method = "pearson")
}




question_1b <- function(datafile){
  
  gender_topic <- chisq.test(datafile$gender, datafile$Topic)
  print(gender_topic)
  gender_performance <- chisq.test(datafile$gender, datafile$Class)
  print(gender_performance)
  #did the code below 3 times and just looked at the general results of the means  
  
  gender_contribution <- split(X <- datafile[ names(datafile) %in% c("raisedhands","Discussion","VisITedResources","gender")], X$gender)
  i = 1
  while(i < 3){
    X <- colMeans(gender_contribution[[i]][sapply(gender_contribution[[i]], is.numeric)], na.rm=TRUE)
    print(gender_contribution[[i]]$gender[1])
    print(X)
    i = i +1
  }
  
  
  #determines mean values for class contribution by results
  numeric_columns <- datafile[ names(datafile) %in% c("raisedhands","AnnouncementsView","Discussion","Class")]
  X <- split(numeric_columns, numeric_columns$Class)
  i =1;
  while(i < 4)
  {
    out <- c(X[[i]]$Class[1],mean(X[[i]]$raisedhands),mean(X[[i]]$AnnouncementsView),mean(X[[i]]$Discussion))
    print(out)
    i = i + 1
  }
    
}

euclidean_dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}
  main()

#girls raise hands more and visit resources more
#
# aggregrate(gender ~ raisedhands, data = datafile, FUN = mean)