#' cLDA Model
#'
#' Obtains a set of filters for labeled time series data
#' so that the between-class distances are maximized, and the within-class
#' distances are minimized.
#'
#' @param Data Matrix of time series on the rows.
#' @param Labels Label of each time series.
#'
#' @return A list containing the filters and their respective importance (g and eig_val), 
#'         the class means (Means), the average of the class means (Mean), and the labels of each class mean (classes). The filters are the columns of the matrix g.
#'
#' @author Grover E. Castro Guzman
#' @author André Fujita
#' 
#' @keywords model
#' 
#' @examples
#' ## Generating 200 time series of length 100 with label 1
#' time_series_signal_1 = sin(matrix(runif(200*100),nrow = 200,ncol = 100))
#' time_series_error_1 = matrix(rnorm(200*100),nrow = 200,ncol = 100)
#' time_series_w_label_1 = time_series_signal_1 + time_series_error_1
#' ## Generating another 200 time series of length 100 with label 2
#' time_series_signal_2 = cos(matrix(runif(200*100),nrow = 200,ncol = 100))
#' time_series_error_2 = matrix(rnorm(200*100),nrow = 200,ncol = 100)
#' time_series_w_label_2 = time_series_signal_2 + time_series_error_2
#' ## Join the time series data in one matrix
#' time_series_data = rbind(time_series_w_label_1,time_series_w_label_2)
#' label_time_series   = c(rep(1,200),rep(2,200))
#' ## obtain the model with the given data
#' clda.model(time_series_data,label_time_series)
#' @export
clda.model <- function(Data,Labels){

  classes <- unique(Labels)
  n_classes   <- rep(0,length(classes))
  Means   <- matrix(0,nrow = length(classes),ncol = ncol(Data))

  for(i in 1:length(classes)){
    index     <- (Labels == classes[i])
    Means[i,] <- colMeans(Data[index,])
    n_classes[i]    <- sum(Labels == classes[i])
  }

  Mean  <- colMeans(Means)
  Dif_Means <- matrix(0,nrow = length(classes),ncol = ncol(Data))

  # Now obtain the between class matrix

  for(i in 1:nrow(Dif_Means)){
    Dif_Means[i,] <- Means[i,] - Mean
  }

  # building the between class matrix
  H_between <- matrix(0,nrow = ncol(Data),ncol = ncol(Data))

  for(i in 1:length(classes)){
    aux <- stats::convolve(Dif_Means[i,],Dif_Means[i,],type = "open")
    lo  <- ncol(Data)
    hi  <- 2*ncol(Data) - 1
    # building the toeplitz matrix
    Toeplitz_aux <- stats::toeplitz(aux[lo:hi])

    H_between <- H_between + (n_classes[i]) * (Toeplitz_aux)
  }

  # Now obtain the within class matrix

  Diff_all <- matrix(0,nrow = nrow(Data),ncol = ncol(Data))

  for(i in 1:nrow(Diff_all)){
    class_index <- which(Labels[i] == classes)
    Diff_all[i,] <- Data[i,] - Means[class_index,]
  }

  # building the between class matrix
  H_within <- matrix(0,nrow = ncol(Data),ncol = ncol(Data))

  for(i in 1:nrow(Diff_all)){
    aux <- stats::convolve(Diff_all[i,],Diff_all[i,],type = "open")
    lo  <- ncol(Data)
    hi  <- 2*ncol(Data) - 1

    # building the toeplitz matrix
    Toeplitz_aux <- stats::toeplitz(aux[lo:hi])

    H_within <- H_within + Toeplitz_aux
  }

  # here the solution ((H_within)^-1) %*% (H_between)

  sol <- MASS::ginv(H_within) %*% H_between

  eigen_sol <- eigen(sol)

  # here the solution
  lo   <- 1
  hi   <- (length(classes) - 1)
  # return the set of filters
  g <- eigen_sol$vectors[,lo:hi]
  # importance of each filter, the first one is the most important
  eig_val  <- eigen_sol$values[lo:hi]

  return (list( "eig_val" = eig_val, "g" = g,
               "Mean" = Mean, "Means" = Means,
               "classes" = classes))
}

#' cLDA classify
#'
#' Classify the time series and obtain the distances between the time series and the centroids of each class.
#'
#' @param model An object returned by the function \code{\link{clda.model}}.
#' @param Data Matrix of time series on the rows.
#'
#'
#' @return A list containing the predicted labels of the time series
#'         and a matrix of distances between the time series and the centroids after applying
#'         the filters obtained by \code{\link{clda.model}}.
#'         
#' @author Grover E. Castro Guzman
#' @author André Fujita
#'
#' @keywords classification
#' @seealso \code{\link{clda.model}}
#' @examples
#'
#' ## Generating 200 time series of length 100 with label 1
#' time_series_signal_1 = sin(matrix(runif(200*100),nrow = 200,ncol = 100))
#' time_series_error_1 = matrix(rnorm(200*100),nrow = 200,ncol = 100)
#' time_series_w_label_1 = time_series_signal_1 + time_series_error_1
#' ## Generating another 200 time series of length 100 with label 2
#' time_series_signal_2 = cos(matrix(runif(200*100),nrow = 200,ncol = 100))
#' time_series_error_2 = matrix(rnorm(200*100),nrow = 200,ncol = 100)
#' time_series_w_label_2 = time_series_signal_2 + time_series_error_2
#' ## Join the time series data in one matrix
#' time_series_data = rbind(time_series_w_label_1,time_series_w_label_2)
#' label_time_series   = c(rep(1,200),rep(2,200))
#' clda_model <- clda.model(time_series_data,label_time_series)
#' ## Create a test set
#' ## data with label 1
#' Data_test_label_1 = sin(matrix(runif(50*100),nrow = 50,ncol = 100))
#' ## data with label 2
#' Data_test_label_2 = cos(matrix(runif(50*100),nrow = 50,ncol = 100))
#' ## join data into a single matrix
#' Data_test = rbind(Data_test_label_1,Data_test_label_2)
#' ## obtain the labels and distances of each time series
#' clda.classify(clda_model,Data_test)
#' @export
clda.classify <- function(model,Data){
  labels   <- rep(-1,nrow(Data))
  distance <- matrix(0,nrow = nrow(Data),ncol = length(model$classes))
  for(i in 1:nrow(Data)){
    # obtain thge label and the distance to each class centroid distance
    time_series <- Data[i,]
    best_class <- -1
    dist  <- Inf
    # distance to each class
    all_dist <- rep(0,length(model$classes))
    for(k in 1:length(model$classes)){
      if(is.null(ncol(model$g))){
        class_mean_filtered  <- stats::convolve(model$Means[k,],model$g,type = "circular")
        time_series_filtered  <- stats::convolve(time_series,model$g,type = "circular")

        # L2 distance
        aux_dist    <- sum((class_mean_filtered - time_series_filtered)^2)
        all_dist[k] <- all_dist[k] + aux_dist
      } else {
        for(j in 1:ncol(model$g)){
          class_mean_filtered  <- stats::convolve(model$Means[k,],model$g[,j],type = "circular")
          time_series_filtered  <- stats::convolve(time_series,model$g[,j],type = "circular")

          # L2 distance
          aux_dist    <- sum((class_mean_filtered - time_series_filtered)^2)
          all_dist[k] <- all_dist[k] + aux_dist
        }
      }
      if(abs(all_dist[k]) < abs(dist)){
        dist  <- all_dist[k]
        best_class <- model$classes[k]
      }
    }
    # saving the best label and the distance to each class centroid
    labels[i] <- best_class
    distance[i,]  <- all_dist
  }

  return (list("labels" = labels,"distance" = distance))
}
