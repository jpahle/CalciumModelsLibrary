library(CalciumModelsLibrary)
context("Basic functionality")

test_that("hello function is present", {
  expect_equal(hello(), "Hello, world!")
  expect_identical(1, 1)
  expect_true(TRUE)
})
