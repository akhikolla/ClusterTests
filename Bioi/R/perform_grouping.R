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
#'
#' @author Zach Colburn
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom stats var
.perform_grouping <- function(
  input,
  critDist,
  use_prog_bar = TRUE
) {
  # Complete type checking was performed at the start of "euclidean_linker"
  # which calls this function. As such, only critical type checks are
  # performed. This is meant mainly to ensure that data passed from the
  # function ".perform_partitioning" which is called by "euclidean_linker" is
  # in the correct format.
  assert_that(class(input)[1] == "matrix")
  assert_that(class(input[1]) %in% c("integer", "numeric"))
  assert_that(nrow(input) >= 1)

  # Return 1 if there is only a single input point.
  if(nrow(input) == 1){
    return(1)
  }

  # Reorder the columns such that the first column accounts for the most
  # variance in position and the last column accounts for the least variance in
  # position.
  column_variances <- apply(input, MARGIN = 2, FUN = var)
  new_column_order <- order(column_variances, decreasing = TRUE)
  input <- input[, new_column_order, drop = FALSE]

  # Order elements in input.
  ordered_elements <- do.call(order, as.data.frame(input))

  # Arrange the rows in input based on the ordering performed above.
  input <- input[ordered_elements,,drop = FALSE]

  # Identify groups.
  groups <- .euclidean_linker_cpp(
    input,
    critDist,
    use_prog_bar = use_prog_bar
  )

  # Change the order of groups in the output to undo the sorting performed
  # above.
  output <- vector(mode = "numeric", length = length(groups))
  output[ordered_elements] <- groups

  # Change group numbers such that they range from 1 to number of groups.
  output <- as.numeric(factor(output))

  # Return the output.
  output
}

