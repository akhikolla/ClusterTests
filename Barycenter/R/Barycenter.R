#'Regularized Wasserstein Barycenters
#'
#' \code{WaBarycenter} takes in a list of matrices representing joint measures on the row and column space and outputs the
#'corresponding Barycenter.
#'The list has to consist of matrices having all the same dimensions, for instance, each matrix represents the normalized weights of the corresponding pixels of images.
#'@author Marcel Klatt
#'@param images A list of matrices satisfying the prerequisites described above.
#'@param maxIter Maximum number of iterations.
#'@param lambda Non-negative regularization parameter (for large lambda the regularized Barycenter is close to its true counterpart). If FALSE the algorithm uses a lambda depending on \code{costm}.
#'@param costm A matrix of pairwise distances between the locations. If FALSE the algorithm uses the usual euclidean distance matrix on a [0,1]x[0,1] equidistant pixel grid.
#'@return The Barycenter of the matrices, represented by a \eqn{n x m} matrix.
#'
#'Given the MNIST dataset, a Barycenter of the digit three is shown below. The Barycenter is based on 4351 images each represented by
#'a 28 x 28 pixel grid, respectively. The values for \code{lambda} and \code{maxIter} were set by default. The dataset is also available in this package (c.f. \link{three}).
#'
#'\figure{threeMNIST.png}{test}
#'@references Cuturi, M.: \code{Fast Computation of Wasserstein Barycenters}, Proceedings of the International Conference on Machine Learning, Beijing, China, 2014
#'@examples #Computation of a Barycenter based on five images representing the digit eight, respectively.
#'WaBarycenter(eight,lambda=10)
#'#For a more reasonable but longer computation!
#'\dontrun{WaBarycenter(eight)}
WaBarycenter <- function(images, maxIter = 10, lambda = FALSE, costm = FALSE){

  time <- proc.time() #to analyze the computation time

  #Check if the specific inputs are correct
    if(is.list(images) == FALSE){
      stop("The images have to be passed as a list each entry representing a matrix!")
    }

    if(length(unique(lapply(images,dim))) == 1){
      dimension <- dim(images[[1]])
    }
    else{
      stop("Dimensions of the images are not equal!")
    }

    if(is.matrix(costm) == FALSE){
      #create a grid of the same dimension as the images on [0,1]² and create the cost matrix
      n <- dimension[1]*dimension[2]
      coord1 <- seq(0,1,length.out = dimension[2])
      coord2 <- rev(seq(0,1,length.out = dimension[1]))
      coordinates <- expand.grid(coord1,coord2)
      costm <- as.matrix(dist(coordinates, diag=TRUE, upper=TRUE))
    }
    else{
      n <- dimension[1]*dimension[2]
      if(identical(dim(costm),rep(n,2)) == FALSE){
      print(costm)
      stop("Dimension of the cost matrix is not compatible with the given images!")
      }
    }
    if(lambda == FALSE){
      lambda <- 60/median(costm)
    }

  ########### Starting the main algorithm ##########

  #initialize a_tild and a_hat
  a_tild <- rep(1/n,n)
  a_hat <- rep(1/n,n)
  t_0 <- 2
  t <- t_0
  #images <- lapply(images,t)


  #iteration using Sinkhorn_Distance (Algorithm 1)
  for(i in 1:maxIter){
    beta <- (t+1)/2
    a <- (1-1/beta) * a_hat + (1/beta) * a_tild
    #Form subgradient with Sinkhorn's Algorithm
    ALPHA <- 0
     for(j in 1:length(images)){
       #This step is based on RcppArmadillo. (Algorithm 3)
       ALPHA <- Subgradient(a,t(images[[j]]),costm,lambda) + ALPHA  #Note, that we need to transpose each matrix, to be compatible with the coordinates defined above.
     }
   # ALPHA <- (1/length(images)) * Reduce("+",lapply(lapply(images,t), Subgradient,a=a,M=costm,lambda))
    ALPHA <- (1/length(images)) * ALPHA
    a_tild <- a_tild*exp(-(t_0)*beta*ALPHA)
    a_tild <- a_tild/sum(a_tild)
    a_hat <- (1-1/beta)*a_hat + (1/beta) * a_tild
    t <- t+1

}

   if(length(unique(dimension)) == 1){
   #Transforming Barycenter s.t. function "image" returns the correct orientation.
   a.temp <- matrix(a,dimension[1],dimension[2],byrow=TRUE)
   a.temp <- a.temp[,nrow(a.temp):1]
   #image(a.temp)
   }

   #Computation time
   print(proc.time()-time)
   #Return the Barycenter
   #result <- matrix(a,dimension[1],dimension[2],byrow=TRUE)
   return(a.temp)
}
