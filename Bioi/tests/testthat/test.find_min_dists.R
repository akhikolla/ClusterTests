context("find_min_dists")

test_that(
  "invalid inputs produce an error",
  {
    expect_error({
      set.seed(10)

      mOne <- as.matrix(data.frame(
        x = rnorm(10),
        y = rbinom(10, 100, 0.5),
        z = runif(10)
      ))

      mTwo <- data.frame(
        x = rnorm(20),
        y = rbinom(20, 100, 0.5),
        z = runif(20)
      )

      find_min_dists(mOne, mTwo)
    })

    expect_error({
      set.seed(10)

      mTwo <- as.matrix(data.frame(
        x = rnorm(10),
        y = rbinom(10, 100, 0.5),
        z = runif(10)
      ))

      mOne <- data.frame(
        x = rnorm(20),
        y = rbinom(20, 100, 0.5),
        z = runif(20)
      )

      find_min_dists(mOne, mTwo)
    })

    expect_error({
      mOne <- TRUE
      mTwo <- "a string"
      euclidean_linker(mOne, mTwo)
    })
  }
)


test_that(
  "outputs are valid",
  {
    expect_equal({
      set.seed(10)

      mOne <- as.matrix(data.frame(
        x = rnorm(10),
        y = rbinom(10, 100, 0.5),
        z = runif(10)
      ))

      mTwo <- as.matrix(data.frame(
        x = rnorm(20),
        y = rbinom(20, 100, 0.5),
        z = runif(20)
      ))

      find_min_dists(mOne, mTwo)
    },
    structure(list(
      dist = c(
        2.10309338977926, 0.17112678621996, 1.07704757389012, 2.36900458966066,
        0.638508938099029, 1.53139116479373, 0.446611676129321,
        0.894671040412461, 1.48945036892257, 0.543050136933722
      ), index = c(11, 12, 9, 4, 12, 15, 16, 10, 13, 15)),
      .Names = c("dist", "index"), row.names = c(NA, -10L), class = "data.frame"
    )
    )
  }
)
