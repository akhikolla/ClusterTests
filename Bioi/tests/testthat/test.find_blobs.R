context("find_blobs")

test_that(
  "invalid inputs produce an error",
  {
    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      find_blobs(mat)
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0.8
      find_blobs(mat[-(1:7),])
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0.8
      find_blobs(mat[,-(1:10)])
    })
  }
)


test_that(
  "outputs are valid",
  {
    expect_equal({
      # Generate a random matrix.
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)

      # Arbitrarily say that everything below 0.8 is background.
      logical_mat <- mat > 0.8

      # Define rownames and column names
      rownames(logical_mat) <- letters[1:7]
      colnames(logical_mat) <- 1:10

      # Find blobs
      find_blobs(logical_mat)
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1, NA, NA,
                NA, NA, NA, NA, 1, NA, NA, 2, 2, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, 3, NA, 4, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, 5, NA, NA, NA, NA, NA), .Dim = c(7L, 10L),
              .Dimnames = list(
                  c("a", "b", "c", "d", "e", "f", "g"),
                  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")))
    )

    expect_equal({
      # Create a logical vector where elements 1, 6-8, and 15-18 are TRUE.
      vec <- rep(FALSE, 20)
      vec[c(1,6:8,15:18)] <- TRUE

      # Add element names.
      names(vec) <- 1:20

      # Find blobs.
      find_blobs(vec)
    },
    structure(c(1, NA, NA, NA, NA, 2, 2, 2, NA, NA, NA, NA, NA, NA,
                3, 3, 3, 3, NA, NA),
              .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                         "11", "12", "13", "14", "15", "16", "17", "18", "19",
                         "20"))
    )


    expect_equal({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0.8
      find_blobs(mat[,1,drop=FALSE])
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA), .Dim = c(7L, 1L))
    )

    expect_equal({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0.8
      find_blobs(mat[1,,drop=FALSE])
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA, 1, NA, NA), .Dim = c(1L, 10L))
    )
  }
)
