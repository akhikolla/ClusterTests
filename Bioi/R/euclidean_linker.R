#' Group PALM/iPALM localizations based on their physical separation distance
#'
#' PALM/iPALM data results in a list of spatial coordinates for fluorophore
#' localizations. This function groups nearby localizations if they are within
#' the provided critical distance from each other.
#'
#' @title Return the group number for each localization.
#'
#' @param input A numeric matrix where each row is a localization and each
#' column is a spatial axis.
#' @param critDist The critical distance for which localizations nearer than
#' this distance are deemed part of the same group.
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
#' @param ... Additional parameters passed to euclidean_linker (i.e.
#' finding_blobs).
#'
#' @author Zach Colburn
#'
#' @examples
#' # Generate random data.
#' #set.seed(10)
#' #input <- as.matrix(data.frame(x=rnorm(10),y=rnorm(10)))
#'
#' # Perform linking.
#' #euclidean_linker(input, 0.4)
#'
#' @export
#'
#' @importFrom assertthat assert_that is.number is.flag
#' @importFrom parallel parLapply makeCluster detectCores
euclidean_linker <- function(
  input,
  critDist,
  use_prog_bar = TRUE,
  run_parallel = FALSE,
  num_cores = NULL,
  partition_req = 5000,
  parallel_call_depth = 3,
  ...
) {
  assert_that(class(input)[1] == "matrix")
  assert_that(class(input[1])[1] %in% c("integer", "numeric"))
  assert_that(nrow(input) >= 1)
  assert_that(nrow(input) < 2147483646)# The C++ max limit for an int minus 1.
  assert_that(ncol(input) > 0)
  assert_that(is.number(critDist))
  assert_that(critDist > 0)
  assert_that(is.flag(use_prog_bar))
  assert_that(is.flag(run_parallel))
  assert_that(
    is.number(partition_req) &&
      (partition_req >= 100) &&
      (as.integer(partition_req) == partition_req)
  )
  assert_that(
    is.null(num_cores) ||
      is.number(num_cores)
  )
  if(run_parallel){use_prog_bar <- FALSE;}
  if(
    is.null(num_cores) && run_parallel
  ){
    num_cores <- detectCores()
  }else if(run_parallel){
    assert_that(
      (as.integer(num_cores) == num_cores) &&
        (num_cores >= 2) &&
        (num_cores <= detectCores())
    )
  }
  finding_blobs <- FALSE
  if("find_blobs" %in% names(match.call())){
    assert_that(is.flag(match.call()[["find_blobs"]]))
    finding_blobs <- match.call()[["find_blobs"]]
  }
  if(finding_blobs){
    min_gap <- 6
  }else{
    min_gap <- 6*critDist
  }

  # Return 1 if there is only a single input point.
  if(nrow(input) == 1){
    return(1)
  }

  # If partitioning:
  if(nrow(input) >= partition_req){
    # This value is used as an identifier of ungrouped points. It is used
    # in ".euclidean_linker_cpp" and ".perform_grouping" as well. The value
    # should not be changed.
    no_group <- -1

    groups <- .perform_partitioning(
      input,
      critDist = critDist,
      use_prog_bar = FALSE,
      run_parallel = run_parallel,
      num_cores = num_cores,
      partition_req = partition_req,
      parallel_call_depth = parallel_call_depth,
      min_gap = min_gap
    )
    output <- as.numeric(factor(groups))

    # Ensure memory from any parallelized threads is retrieved.
    gc()

    return(output)
  }

  # If not partitioning:
  output <- .perform_grouping(
    input,
    critDist,
    use_prog_bar = use_prog_bar
  )

  # Return the output.
  output
}

