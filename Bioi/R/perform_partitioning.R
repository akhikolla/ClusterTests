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
#' @param min_gap The minimum width of any dimension created during
#' partitioning.
#'
#' @author Zach Colburn
#'
#' @export
#'
#' @importFrom assertthat assert_that is.number is.flag
#' @importFrom parallel parLapply makeCluster detectCores
#' @importFrom igraph components graph.data.frame
#' @importFrom dplyr mutate_all
.perform_partitioning <- function(
  input,
  critDist,
  use_prog_bar = TRUE,
  run_parallel = FALSE,
  num_cores = NULL,
  partition_req = 5000,
  parallel_call_depth = 3,
  min_gap = NULL
) {
  # Perform essential type checking.
  assert_that((class(input)[1] == "matrix"))
  assert_that(is.number(min_gap))

  # Define the ungrouped value. Don't change this.
  no_group <- -1

  # Check whether parallel_call_depth is 1. If so set run_parallel to FALSE.
  if(parallel_call_depth == 1){
    run_parallel <- FALSE
  }

  # Handle the case of only a single input point.
  if(nrow(input) == 1){
    return(1)
  }

  # Specify partitions.
  use_prog_bar <- FALSE
  ranges <- apply(input, 2, function(item){max(item) - min(item)})
  dim_to_split <- which.max(ranges)

  # If the partition is too small, then run euclidean linker on this data set
  # without performing partitioning, regardless of how many points are in the
  # partition. Also, if the number of points in this partition is less than the
  # partition requirement, then perform grouping.
  if(
    (min(ranges) <= min_gap) ||
    (nrow(input) <= partition_req)
  ){
    groups <- .perform_grouping(
      input,
      critDist,
      use_prog_bar = use_prog_bar
    )
    return(groups)
  }

  # Determine partitioning bounds.
  max_pos <- max(input[,dim_to_split])
  min_pos <- min(input[,dim_to_split])
  mid_pos <- min_pos + (max_pos - min_pos)/2

  # Assign points to partitions and stem group.
  partition <- rep(1, nrow(input))
  partition[input[,dim_to_split] > mid_pos] <- 2
  stem <- (input[,dim_to_split] >= (mid_pos - min_gap/2)) &
    (input[,dim_to_split] < (mid_pos + min_gap/2))

  # If the stem and both partitions don't all have points then group now.
  if(!(
    (sum(stem) >= 1) &&
    (sum(partition == 1) >= 1) &&
    (sum(partition == 2) >= 1)
  )){
    groups <- .perform_grouping(
      input,
      critDist,
      use_prog_bar = use_prog_bar
    )
    groups <- as.numeric(factor(groups))
    return(groups)
  }

  # Perform grouping.
  if(!run_parallel){
    groups <- lapply(
      list(
        input[partition == 1,, drop = FALSE],
        input[partition == 2,, drop = FALSE],
        input[stem,, drop = FALSE]
      ),
      function(item){
        .perform_partitioning(
          item,
          critDist,
          use_prog_bar = FALSE,
          run_parallel = run_parallel,
          num_cores = NULL,
          partition_req = partition_req,
          min_gap = min_gap
        )
      }
    )
  }else{
    # Perform operations in parallel.
    groups <- parLapply(
      makeCluster(num_cores),
      list(
        input[partition == 1,, drop = FALSE],
        input[partition == 2,, drop = FALSE],
        input[stem,, drop = FALSE]
      ),
      function(
        item,
        critDist,
        run_parallel,
        num_cores,
        partition_req,
        parallel_call_depth,
        .perform_partitioning,
        min_gap
      ){
        .perform_partitioning(
          item,
          critDist,
          use_prog_bar = FALSE,
          run_parallel = run_parallel,
          num_cores = num_cores,
          partition_req = partition_req,
          parallel_call_depth = parallel_call_depth - 1,
          min_gap = min_gap
        )
      }, critDist,
      run_parallel,
      num_cores,
      partition_req,
      parallel_call_depth,
      .perform_partitioning,
      min_gap
    )
  }

  # Shift up par 2 group numbers so they don't overlap with par 1 group
  # numbers.
  groups[[2]] <- groups[[2]] + max(groups[[1]])

  # Create a matrix of these values.
  mat <- matrix(no_group, nrow = length(partition), ncol = 3)
  mat[partition == 1, 1] <- groups[[1]]
  mat[partition == 2, 2] <- groups[[2]]
  mat[stem, 3] <- groups[[3]]

  # Set the values in the stem column to be larger than the values in the other
  # columns.
  mat[mat[,3] != -1, 3] <- mat[mat[,3] != -1, 3] + max(mat[, c(1,2)])

  # If there are no points in the stem, then retrieve the final group numbers.
  if(length(unique(mat[mat[,3] != no_group, 3])) == 0){
    mat <- mat[, 1:2, drop = FALSE]
    groups <- t(mat)
    groups <- groups[groups != no_group]
    groups <- as.numeric(factor(groups))
    return(groups)
  }

  # Structure incoming data.
  colnames(mat) <- c("pOne", "pTwo", "stem")
  mat <- as.data.frame(mat)

  # Link pOne to stem and update pTwo accordingly.
  ## Select pOne and stem.
  pOne <- mat[,c("pOne", "stem"), drop = FALSE]
  ## Determine the largest group number.
  maxGroup <- max(apply(mat, 2, max))
  ## Create a data.frame of the links to be made.
  pOne_dict <- pOne[
    (pOne[,"pOne"] != no_group) &
      (pOne[,"stem"] != no_group),
    , drop = FALSE
    ]
  ## Convert those values to factors.
  pOne_dict <- mutate_all(pOne_dict, as.factor)
  ## Perform linking.
  members <- components(
    graph.data.frame(
      pOne_dict,
      directed = FALSE
    ),
    mode = "weak"
  )$membership
  ## Convert mat to a matrix in order then update all group numbers accordingly.
  mat <- as.matrix(mat)
  logicalVector <- (mat != no_group) & (as.character(mat) %in% names(members))
  mat[logicalVector] <- members[as.character(mat[logicalVector])] + maxGroup
  mat <- as.data.frame(mat)

  # Link pOne to stem and update pTwo accordingly.
  ## Select pOne and stem.
  pTwo <- mat[,c("pTwo", "stem"), drop = FALSE]
  ## Determine the largest group number.
  maxGroup <- max(apply(mat, 2, max))
  ## Create a data.frame of the links to be made.
  pTwo_dict <- pTwo[
    (pTwo[,"pTwo"] != no_group) &
      (pTwo[,"stem"] != no_group),
    , drop = FALSE
    ]
  ## Convert those values to factors.
  pTwo_dict <- mutate_all(pTwo_dict, as.factor)
  ## Perform linking.
  members <- components(
    graph.data.frame(
      pTwo_dict,
      directed = FALSE
    ),
    mode = "weak"
  )$membership
  ## Convert mat to a matrix then update all group numbers accordingly.
  mat <- as.matrix(mat)
  logicalVector <- (mat != no_group) & (as.character(mat) %in% names(members))
  mat[logicalVector] <- members[as.character(mat[logicalVector])] + maxGroup

  # Select the pOne and pTwo columns
  mat <- mat[,c("pOne", "pTwo"), drop = FALSE]
  # Transpose the matrix and then remove all no_group values.
  mat <- t(mat)
  groups <- mat[mat != no_group]

  gc()
  as.numeric(factor(groups))
}

