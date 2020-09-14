#' @title For all points in matrix 1, return the distance to and index of the
#' nearest point in matrix 2.
#'
#' @description
#' Find the shortest distance between each point in one data set and the points
#' in a second set.
#'
#' This function determines the distance between every point in mOne and the
#' nearest point in mTwo.
#'
#' @param mOne A numeric matrix where each row is a localization and each
#' column is a spatial axis.
#' @param mTwo A numeric matrix with the same number of columns as mOne.
#'
#' @author Zach Colburn
#'
#' @examples
#' # Generate random data.
#' set.seed(10)
#'
#' mOne <- as.matrix(data.frame(
#' x = rnorm(10),
#' y = rbinom(10, 100, 0.5),
#' z = runif(10)
#' ))
#'
#' mTwo <- as.matrix(data.frame(
#' x = rnorm(20),
#' y = rbinom(20, 100, 0.5),
#' z = runif(20)
#' ))
#'
#' # Find the minimum distance between each point in mOne and the points in
#' # mTwo.
#' find_min_dists(mOne, mTwo)
#'
#' @export
#'
#' @importFrom assertthat assert_that
find_min_dists <- function(mOne, mTwo) {
  # Perform type checking.
  assert_that(class(mOne)[1] == "matrix")
  assert_that(length(mOne) >= 1)
  assert_that(class(mOne[1]) %in% c("integer", "numeric"))
  assert_that(nrow(mOne) >= 1)
  assert_that(ncol(mOne) > 0)

  assert_that(class(mTwo)[1] == "matrix")
  assert_that(length(mTwo) >= 1)
  assert_that(class(mTwo[1]) %in% c("integer", "numeric"))
  assert_that(nrow(mTwo) >= 1)
  assert_that(ncol(mTwo) > 0)

  # Find the minimum distance between each point in mOne and the points in mTwo.
  .find_min_dists_cpp(mOne, mTwo)
}
