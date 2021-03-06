################################################################################

context("RANDOM_SVD")

set.seed(SEED)

################################################################################

TOL <- 1e-4

# function for sampling scaling
sampleScale <- function() {
  tmp <- sample(list(c(TRUE, FALSE),
                     c(TRUE, TRUE),
                     c(FALSE, FALSE)))[[1]]
  list(center = tmp[1], scale = tmp[2])
}

# Simulating some data
N <- 73
M <- 43
x <- matrix(rnorm(N * M, mean = 100, sd = 5), N)

################################################################################

test_that("equality with prcomp", {
  for (t in TEST.TYPES) {
    X <- `if`(t == "raw", asFBMcode(x), big_copy(x, type = t))

    k <- sample(c(2, 5, 20), 1) # 2, 5 or 20

    test <- big_randomSVD(X, k = k, tol = 1e-10, ncores = test_cores())
    pca <- prcomp(X[], center = FALSE, scale. = FALSE)
    expect_equal(diffPCs(predict(test), pca$x), 0, tolerance = TOL)
    expect_equal(diffPCs(test$v, pca$rotation), 0, tolerance = TOL)

    sc <- sampleScale()
    test <- big_randomSVD(X, k = k, tol = 1e-10, ncores = test_cores(),
                          fun.scaling = big_scale(center = sc$center,
                                                  scale = sc$scale))
    pca <- prcomp(X[], center = sc$center, scale. = sc$scale)
    expect_equal(diffPCs(predict(test), pca$x), 0, tolerance = TOL)
    expect_equal(diffPCs(test$v, pca$rotation), 0, tolerance = TOL)
    if (sc$center) expect_equal(test$center, pca$center)
    if (sc$scale)  expect_equal(test$scale,  pca$scale)

    p <- plot(test, type = sample(c("screeplot", "scores", "loadings"), 1))
    expect_s3_class(p, "ggplot")

    expect_error(predict(test, abc = 2), "Argument 'abc' not used.")
    expect_error(plot(test, abc = 2), "Argument 'abc' not used.")
  }
})

################################################################################

test_that("equality with prcomp with half of the data", {
  ind <- sample(N, N / 2)
  ind2 <- setdiff(1:N, ind)

  for (t in TEST.TYPES) {
    X <- `if`(t == "raw", asFBMcode(x), big_copy(x, type = t))

    k <- sample(c(2, 5, 20), 1) # 2, 5 or 20
    sc <- sampleScale()

    test <- big_randomSVD(X, ncores = test_cores(),
                          ind.row = ind, k = k, tol = 1e-10,
                          fun.scaling = big_scale(center = sc$center,
                                                  scale = sc$scale))
    pca <- prcomp(X[ind, ], center = sc$center, scale. = sc$scale)

    expect_equal(diffPCs(predict(test), pca$x), 0, tolerance = TOL)
    expect_equal(diffPCs(test$v, pca$rotation), 0, tolerance = TOL)

    if (sc$center) expect_equal(test$center, pca$center)
    if (sc$scale)  expect_equal(test$scale,  pca$scale)

    expect_equal(diffPCs(predict(test, X, ind.row = ind2),
                         predict(pca, X[ind2, ])), 0, tolerance = TOL)

    p <- plot(test, type = sample(c("screeplot", "scores", "loadings"), 1))
    expect_s3_class(p, "ggplot")
  }
})

################################################################################

test_that("equality with prcomp with half of half of the data", {
  ind <- sample(N, N / 2)
  ind2 <- setdiff(1:N, ind)
  ind.col <- sample(M, M / 2)
  for (t in TEST.TYPES) {
    X <- `if`(t == "raw", asFBMcode(x), big_copy(x, type = t))

    k <- sample(c(2, 5, 20), 1) # 2, 5 or 20
    sc <- sampleScale()

    test <- big_randomSVD(X, ind.row = ind, ind.col = ind.col,
                          k = k, tol = 1e-10, ncores = test_cores(),
                          fun.scaling = big_scale(center = sc$center,
                                                  scale = sc$scale))
    pca <- prcomp(X[ind, ind.col], center = sc$center, scale. = sc$scale)

    expect_equal(diffPCs(predict(test), pca$x), 0, tolerance = TOL)
    expect_equal(diffPCs(test$v, pca$rotation), 0, tolerance = TOL)

    if (sc$center) expect_equal(test$center, pca$center)
    if (sc$scale)  expect_equal(test$scale,  pca$scale)

    expect_equal(diffPCs(predict(test, X, ind.row = ind2, ind.col = ind.col),
                         predict(pca, X[ind2, ind.col])), 0, tolerance = TOL)

    p <- plot(test, type = sample(c("screeplot", "scores", "loadings"), 1))
    expect_s3_class(p, "ggplot")
  }
})

################################################################################
