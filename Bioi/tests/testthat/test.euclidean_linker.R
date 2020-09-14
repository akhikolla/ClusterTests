context("euclidean_linker")

test_that(
  "invalid inputs produce an error",
  {
    expect_error({
      input <- data.frame(x = 1:10, y = 1:10)
      crit_dist <- 1.5
      euclidean_linker(input, crit_dist)
    })

    expect_error({
      input <- as.matrix(data.frame(x = 1:10, y = 1:10))
      crit_dist <- "TRUE"
      euclidean_linker(input, crit_dist)
    })

    expect_error({
      input <- as.matrix(data.frame(x = 1:10, y = 1:10))
      crit_dist <- TRUE
      euclidean_linker(input, crit_dist)
    })
  }
)


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

    expect_true({
      set.seed(10)
      input <- as.matrix(data.frame(x = rnorm(10), y = rnorm(10)))
      crit_dist <- 0.4
      result <- euclidean_linker(input, crit_dist)
      (result[2] == tail(result, 1) && (length(unique(result))  == 9))
    })

    expect_true({
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
      length(unique(result)) == 3
    })

    expect_equal({
      mat <- matrix(TRUE, nrow = 501, ncol = 501)
      result <- Bioi::find_blobs(mat)
      unique(as.vector(result))
    },
    1
    )
  }
)


