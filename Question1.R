# install.packages('igraph')
# install.packages('cccd')
# install.packages('caret')
# install.packages('Boruta')
main <- function() {
  library('cccd')
  library('igraph')
  library('Boruta')
  # datafile <- read.csv('../Datasets/xAPI-edu-Data.csv')
  # pearsons <- question_1a(datafile)
  # question_1b(datafile)
  # datafile <- read.csv("../Datasets/AlzheimersDisease.csv", header = TRUE)
  datafile2 <- read.csv("../Datasets/AlzheimersDisease.csv", header = FALSE)
  
  # question_2(datafile)
  question_4()
  # question_3(datafile2)
  # generate_sheetAlzheimersCSV()
  # # generate_alzheimers_matrixes(datafile)
  # get_relative_neighbourhoods()
  # question_3(datafile) 
}

question_4 <- function(){
  datafile <- read.csv("../Datasets/Iris.csv", header = TRUE)
  rownames(datafile) <- datafile[,1]
  datafile[,1] <- NULL
  checking <- datafile[1]
  class <-  datafile[5]
  datafile[,5] <- NULL
  # generate_distances(datafile,"IrisRowsMatrix","IrisAttributesMatrix")
  datafile[,1] <- class
  datafile[,5] <- checking
  colnames(datafile)[1] <- "CLASS"
  colnames(datafile)[5] <- "SepalLengthCm"
  print(datafile)
  feature_selection(datafile,"IrisFeatures",FALSE)
}
question_2 <- function(datafile){

  rownames(datafile) <- datafile[,1]
  datafile[,1] <- NULL
  print(datafile)
  generate_distances(datafile,"samplesMatrix","ProteinsMatrix")
}
question_3 <- function(datafile){
  datafile <- add_rownames(datafile)
  datafile <- t(datafile)
  # print(datafile)
  #prints classes that are important in determining Alzheimers
  proteins <- feature_selection(datafile,"proteins",TRUE)
  #uses those printed classes to do the classification
  # generate_sheetAlzheimersCSV()
}

add_rownames <- function(datafile){
  rownames(datafile) <- datafile[,1]
  datafile[,1] <- NULL
  rownames(datafile)[1]<- "CLASS"
  datafile
}

feature_selection <- function(data,filename, alzheimers){
  set.seed(111)
  bor <- Boruta(as.factor(data[,1]), x = data)
  print(bor)
  bor <- TentativeRoughFix(bor, averageOver = Inf)
  x <- getSelectedAttributes(bor, withTentative = FALSE)
  x <- as.data.frame(x)
  print(x)
  a = 1
  if(alzheimers == TRUE)
  {
    while(a <= nrow(x))
    {
      x[a,1] <- gsub(".","-",x[a,1], fixed = TRUE)
      a = a + 1
    }
  }
  
  matrix1 <- matrix(data = data[,colnames(data) %in% x[,1]], nrow = nrow(data), ncol = nrow(x))
  colnames(matrix1) <- x[,1]
  rownames(matrix1) <- rownames(data)
  a = 1
  if(alzheimers == TRUE)
  {
    matrix1 <- fix_proteins(matrix1)
  }
  name = paste("Results/",filename,".csv", sep = "")
  write.csv(matrix1,file = name)
  x
}

fix_proteins <-function(matrix){
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


#creates RNG graphs with a matrix and a title given
relative_neighbourhoods <- function(matrix,title){
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

generate_distances <- function(datafile, filename1,filename2){
  row_matrix <- get_distance_rows(datafile)
  cols_matrix <- get_distance_cols(datafile)
  name1 = paste("Results/",filename1,"distancematrix.csv", sep = "")
  name2 = paste("Results/",filename2,"distancematrix.csv", sep = "")
  write.csv(row_matrix,file = name1)
  write.csv(cols_matrix,file = name2)
  relative_neighbourhoods(row_matrix,filename1)
  relative_neighbourhoods(cols_matrix,filename2)
  
  rows_graph <- mst(graph_from_adjacency_matrix(row_matrix, mode = "undirected",weighted = TRUE))
  cols_graph <- mst(graph_from_adjacency_matrix(cols_matrix, mode = "undirected",weighted = TRUE))
  
  V(rows_graph)$label <- colnames(row_matrix)
  E(rows_graph)$label <- round(E(rows_graph)$weight, 3)
  V(cols_graph)$label <- colnames(cols_matrix)
  E(cols_graph)$label <- round(E(cols_graph)$weight, 3)
  print(sum(E(rows_graph)$weight))
  print(sum(E(cols_graph)$weight))
  name1 = paste("Results/",filename1,"MST.gml", sep = "")
  name2 = paste("Results/",filename2,"MST.gml", sep = "")
   
  write_graph(rows_graph, name1,format = 'gml')
  write_graph(cols_graph,name2,format = 'gml')
}

#generates euclidean distance matrix of rows
get_distance_rows <- function(datafile)
{
  rows <- matrix(nrow = nrow(datafile), ncol = nrow(datafile))
  colnames(rows) <- rownames(datafile)
  rownames(rows) <- rownames(datafile)
  i <- 1
  a <- 1
  print(ncol(rows))
  while(i < nrow(rows)+1){
    while(a < ncol(rows)+1)
    {
      rows[i,a] <- euclidean_dist(datafile[i,],datafile[a,])
      a = a + 1
    }
    i = i + 1
    a = 1
    print("loop done")
  }
  print("rows done")
  return(rows)
}

#generates euclidean distance matrix of columns
get_distance_cols <- function(datafile)
{
  
  columns <- matrix(nrow = ncol(datafile), ncol = ncol(datafile))
  colnames(columns) <- colnames(datafile)
  rownames(columns) <- colnames(datafile)
  print(datafile)
  i <- 1
  a <- 1
  while(i < nrow(columns)+1){
    while(a < ncol(columns)+1)
    {
      columns[i,a] <- euclidean_dist(datafile[,i],datafile[,a])
      a = a + 1
    }
    i = i + 1
    a = 1
  }
  print("distance col done")
  return(columns)
}

generate_sheetAlzheimersCSV <- function(){
  test_ad <- read.csv("../Datasets/TestSetAD.csv", header = FALSE)
  test_ad <- add_rownames(test_ad)
  test_ad <- t(test_ad)
  test_ad <- fix_proteins(test_ad)
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
  test_mci <- fix_proteins(test_mci)
  write.csv(test_mci,"testMCI.csv")
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

# aggregrate(gender ~ raisedhands, data = datafile, FUN = mean)