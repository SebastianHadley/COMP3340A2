# install.packages('igraph')
# install.packages('cccd')
# install.packages('caret')
# install.packages('Boruta')
main <- function() {
  library('cccd')
  library('igraph')
  library('Boruta')
  datafile <- read.csv('../Datasets/xAPI-edu-Data.csv')
  pearsons <- question_1a(datafile)
  question_1b(datafile)
  datafile <- read.csv("../Datasets/AlzheimersDisease.csv", header = FALSE)
  # question_3(datafile)
  generate_sheetAlzheimersCSV()
  # # generate_matrixes(datafile)
  # get_relative_neighbourhoods()
  # question_3(datafile) 
}
 
generate_sheetAlzheimersCSV <- function(){
  test_ad <- read.csv("../Datasets/TestSetAD.csv", header = FALSE)
  test_ad <- add_rownames(test_ad)
  test_ad <- t(test_ad)
  test_ad <- fix_classes(test_ad)
  write.csv(test_ad,"testAD.csv")
  test_mci <- read.csv("../Datasets/TestSetMCI.csv", header = TRUE)
  colnames(test_mci) <- NULL
  temp <- test_mci[1,]
  temp[1] <- "CLASS"
  test_mci[1,] <- temp
  x <- nrow(test_mci) - 121
  test_mci <- head(test_mci, -x)
  test_mci <- t(test_mci)
  colnames(test_mci) <-  test_mci[1,]
  test_mci <- test_mci[-c(1),]
  test_mci <- fix_classes(test_mci)
  write.csv(test_mci,"testMCI.csv")
}

question_3 <- function(datafile){
  datafile <- add_rownames(datafile)
  datafile <- t(datafile)
  #prints classes that are important in determining Alzheimers
  proteins <- feature_selection(datafile)
  #uses those printed classes to do the classification
  # classification_code(proteins)
}

add_rownames <- function(datafile)
{
  rownames(datafile) <- datafile[,1]
  datafile[,1] <- NULL
  rownames(datafile)[1]<- "CLASS"
  datafile
}

feature_selection <- function(data){
  set.seed(111)
  bor <- Boruta(as.factor(data[,1]), x = data)
  bor <- TentativeRoughFix(bor, averageOver = Inf)
  x <- getSelectedAttributes(bor, withTentative = FALSE)
  x <- as.data.frame(x)
  a = 1
  while(a <= nrow(x))
  {
    x[a,1] <- gsub(".","-",x[a,1], fixed = TRUE)
    a = a + 1
  }
  matrix1 <- matrix(data = data[,colnames(data) %in% x[,1]], nrow = nrow(data), ncol = nrow(x))
  colnames(matrix1) <- x[,1]
  rownames(matrix1) <- rownames(data)
  a = 1
  matrix1 <- fix_classes(matrix1)
  write.csv(matrix1, file = "Results/features.csv")
  x
}

fix_classes <-function(matrix){
  a = 1
  while(a <= nrow(matrix))
  {
    if(matrix[a,1] != 'AD'){
      matrix[a,1] = "NON_AD"
    }
    a = a + 1
  }
  matrix
}
get_relative_neighbourhoods <- function(){
  samples_matrix <- read.csv("../Datasets/samples_distance_matrix.csv", header = TRUE)
  proteins_matrix <- read.csv("../Datasets/proteins_distance_matrix.csv",header= TRUE)
  rownames(samples_matrix) <- samples_matrix[,1]
  rownames(proteins_matrix) <- proteins_matrix[,1]
  proteins_matrix[,1] <- NULL
  samples_matrix[,1] <- NULL
  relative_neighbourhoods(samples_matrix,"samples")
  relative_neighbourhoods(proteins_matrix,"proteins")
}


relative_neighbourhoods <- function(matrix,title)
{
  matrix1 <- as.matrix(matrix)
  rng_graph <- rng(dx = matrix1,r = 1, algorithm = 'cover_tree')
  rng_graph <- as.undirected(rng_graph)
  #adds labels
  V(rng_graph)$label <- colnames(matrix1)
  E(rng_graph)$weight <- apply(get.edges(rng_graph, 1:gsize(rng_graph)),1,function(x)matrix1[x[1],x[2]])
  E(rng_graph)$label <- round(E(rng_graph)$weight, 3)
  name <- paste("Results/",title,"RNG.gml", sep ="")
  plot(rng_graph)
  write_graph(rng_graph, name, format = 'gml')
}


generate_matrixes <- function(datafile)
{
  rownames(datafile) <- datafile[,1]
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