context("euclidean_linker_cpp")


test_that(
  "outputs are valid",
  {
    expect_true({
      vec_length <- 10000
      input <- as.matrix(data.frame(x = 1:vec_length, y = 1:vec_length))
      crit_dist <- 1.5
      result <- euclidean_linker(input, crit_dist)
      all(result == 1)
    })

    expect_equal({
      set.seed(10)
      input <- as.matrix(data.frame(x = rnorm(10), y = rnorm(10)))
      crit_dist <- 0.4
      euclidean_linker(input, crit_dist)
    },
    c(7, 6, 2, 4, 8, 9, 3, 5, 1, 6)
    )

    expect_equal({
      set.seed(10)
      num <- 1000
      input <- as.matrix(data.frame(
        x = runif(num)*40,
        y = runif(num)*40,
        z = runif(num)*40
      ))
      crit_dist <- 5
      result <- euclidean_linker(
        input,
        crit_dist,
        partition_req = 200,
        run_parallel = FALSE
      )
      unique(result)
    },
    c(3, 2, 1)
    )

    expect_equal({
      mat <- matrix(TRUE, nrow = 501, ncol = 501)
      result <- Bioi::find_blobs(mat)
      unique(as.vector(result))
    },
    1
    )
  }
)