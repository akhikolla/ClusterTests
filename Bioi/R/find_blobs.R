#' @title Assign all neighboring pixels the same group number.
#'
#' @description
#' Perform connected-component labeling to group continuous, thresholded
#' objects in 3-dimensional arrays.
#'
#' This function takes a vector, matrix, or 3-dimensional array where each
#' element is TRUE if it corresponds to an object-positive index or FALSE if it
#' corresponds to a background index. An object of the same dimension as the
#' input is returned. All connected object indices take the value of their
#' group number and all background indices take the value NA.
#'
#' @param arr A vector, matrix, or 3-dimensional array where object-positive
#' elements are denoted by the value TRUE and background elements are denoted
#' by the value FALSE.
#' @param use_prog_bar TRUE/FALSE indicating whether a progress bar should be
#' used. This is only available when run_parallel is FALSE.
#' @param run_parallel TRUE/FALSE indicating whether operations should be
#' performed in parallel. This is only valid if partitioning is performed.
#' @param num_cores The number of cores to use if running in parallel.
#' @param partition_req The minimum number of points required to create a new
#' partition.
#' @param parallel_call_depth The number of levels of partitioning that should
#' be performed before terminating calls to run operations in parallel. The
#' number of threads opened when running in parallel is equal to
#' 2^(parallel_call_depth)*num_cores.
#'
#' @author Zach Colburn
#'
#' @examples
#' # Generate a random matrix.
#' set.seed(10)
#' mat <- matrix(runif(70), nrow = 7)
#'
#' # Arbitrarily say that everything below 0.8 is background.
#' logical_mat <- mat > 0.8
#'
#' # Find blobs.
#' find_blobs(logical_mat)
#'
#' @export
#'
#' @importFrom assertthat assert_that noNA is.number
find_blobs <- function(
  arr,
  use_prog_bar = TRUE,
  run_parallel = FALSE,
  num_cores = NULL,
  partition_req = NULL,
  parallel_call_depth = 3
) {
  # Perform type checking.
  assert_that(is.vector(arr) || is.matrix(arr) || is.array(arr))
  assert_that(length(arr) >= 1)
  assert_that(class(arr[1])[1] == "logical")
  assert_that(noNA(arr))
  if(is.null(partition_req)){partition_req <- 1000000}
  assert_that(
    is.number(partition_req) &&
      (as.integer(partition_req) == partition_req) &&
      (partition_req > 0)
  )

  # Get object class and attributes.
  initial_class <- class(arr)[1]
  arr_attributes <- attributes(arr)

  # Convert arr to an array and store object indices in "input".
  arr <- as.array(arr)
  input <- as.matrix(which(arr, arr.ind = TRUE))

  # Initialize the output object.
  output <- array(NA, dim = dim(arr), dimnames = dimnames(arr))

  # If there are no object indices then return output without performing
  # any more operations.
  if(nrow(input) == 0){
    print("There are no objects to link!")
    return(output)
  }

  # Use the euclidean_linker_cpp function to link neighboring object indices.
  # Since the critical distance is sqrt(3), every neighboring object index
  # (both horizontally/vertically and diagonally) will be joined in 1, 2, or 3
  # dimensions. "links" is a vector of group numbers.
  links <- euclidean_linker(
    input,
    critDist = sqrt(3),
    use_prog_bar = use_prog_bar,
    run_parallel = run_parallel,
    num_cores = num_cores,
    partition_req = partition_req,
    parallel_call_depth = parallel_call_depth,
    find_blobs = TRUE
  )

  # Assign indices in the output object their respective group number.
  output[input] <- links

  # Convert the output object to the class of the original arr object.
  if(initial_class == "logical"){
    output <- as.vector(output)
  }else if(initial_class == "matrix"){
    output <- as.matrix(output)
  }

  # Restore the object's attributes. This is mainly for restoring vector
  # element names and matrix row/column names).
  attributes(output) <- arr_attributes

  # Return the output object.
  return(output)
}


#' @title Assign all neighboring pixels the same group number.
#'
#' @description
#' This function is deprecated. It now calls the more efficient find_blobs
#' method.
#'
#' This function takes a matrix corresponding to a thresholded image and
#' returns a matrix of the same size, where all adjacent, thresholded pixels
#' are the same integer corresponding to that object's cluster number.
#'
#' @param img A thresholded matrix (where non-object pixels are assigned a
#' value of 0).
#' @param pixRange This parameter is now obsolete. Previously, the parameter
#' denoted an integer number of pixels to specify a search region. Execution
#' was faster when this value was small. However, the value needed to be larger
#' than the diameter of the largest continuous object in the image.
#'
#' @author Zach Colburn
#'
#' @examples
#' # Generate a random matrix.
#' set.seed(10)
#' mat <- matrix(runif(70), nrow = 7)
#'
#' # Arbitrarily say that everything below 0.8 is background.
#' mat[mat < 0.8] <- 0
#'
#' # Find blobs.
#' identify_thresholded_objects(mat)
#'
#' @export
#'
#' @importFrom assertthat assert_that is.number
identify_thresholded_objects <- function(img, pixRange = 50){
  # Perform type checking. This function was meant to receive different inputs
  # than the function find_blobs which replaces it. Type checking is performed
  # to ensure object inputs are backwards compatible.
  assert_that(class(img)[1] == "matrix")
  assert_that(length(img) >= 1)
  assert_that(class(img[1]) %in% c("integer", "numeric"))
  assert_that(nrow(img) >= 1)
  assert_that(ncol(img) > 0)

  # If pixRange has changed, then inform the user that its use is deprecated.
  assert_that(is.number(pixRange))
  if(pixRange != 50){
    print("The use of pixRange is deprecated.")
  }

  # Convert the input to a logical matrix.
  img <- img != 0

  # Perform connected component labeling.
  find_blobs(img)
}