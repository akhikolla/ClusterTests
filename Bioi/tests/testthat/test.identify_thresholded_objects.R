context("identify_thresholded_objects")

test_that(
  "invalid inputs produce an error",
  {
    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0
      identify_thresholded_objects(mat)
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat[-(1:7),])
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat[,-(1:10)])
    })
  }
)

test_that(
  "outputs are valid",
  {
    expect_equal({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat)
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1, NA, NA,
                NA, NA, NA, NA, 1, NA, NA, 2, 2, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, 3, NA, 4, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, 5, NA, NA, NA, NA, NA), .Dim = c(7L, 10L))
    )
  }
)


test_that(
  "use of pixRange causes no problems",
  {
    expect_equal({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat, 51)# 51 beacuse pixRange defaults to 50.
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1, NA, NA,
                NA, NA, NA, NA, 1, NA, NA, 2, 2, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, 3, NA, 4, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, 5, NA, NA, NA, NA, NA), .Dim = c(7L, 10L))
    )
  }
)