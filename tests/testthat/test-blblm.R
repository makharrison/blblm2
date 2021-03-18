test_that("output of blblm is of class blblm", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
  expect_identical(class(fit), "blblm")
})

test_that("coefficients can be found", {
  blb <- coef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE))
  notblb <- lm(mpg ~ wt * hp, data = mtcars)$coefficients
  expect_identical(length(blb), length(notblb))
})

test_that("coefficients are numeric", {
  blb <- coef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE))
  expect_identical(class(blb), "numeric")
})

test_that("sigma works", {
  blb <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
  expect_identical(class(sigma(blb)), "numeric")
})

test_that("parallel works", {
  library(furrr)
  plan(multiprocess, workers = 2)
  blb1 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)
  blb2 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
  expect_identical(coef(blb1), coef(blb2), tolerance = 10)
})

test_that("class of blblog object", {
  l <- blblog(vs ~ wt + disp, data = mtcars, m = 3, B = 100)
  expect_identical(class(l), "blblog")
})
