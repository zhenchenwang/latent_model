# Unit tests for the functions in latentModel.R

# Load the testthat library
library(testthat)

# Source the script with the functions to be tested
# Using a relative path from the project root.
# This assumes the tests are run from the project root directory.
source("scripts/latentModel.R")

# The tests will be added in the next steps.


# -----------------------------------------------------------------------------
# Tests for Helper Functions
# -----------------------------------------------------------------------------

# Create a small sample dataframe for testing
test_data <- data.frame(
  age = c(50, 60, 70),
  bmi = c(25.1, 30.2, 22.3),
  strokeha = c(0, 1, 0),
  gender = c("M", "F", "M")
)

test_that("format_variables correctly converts column types", {
  formatted_data <- format_variables(test_data)

  # Check if numeric columns are numeric
  expect_true(is.numeric(formatted_data$age))
  expect_true(is.numeric(formatted_data$bmi))

  # Check if factor columns are factors
  expect_true(is.factor(formatted_data$strokeha))
  expect_true(is.factor(formatted_data$gender))
})

test_that("discretize_data correctly discretizes numeric columns", {
  discretization_bins <- list(age = 2, bmi = 3)
  discretized_data <- discretize_data(test_data, discretization_bins)

  # Check if the specified columns are now factors
  expect_true(is.factor(discretized_data$age))
  expect_true(is.factor(discretized_data$bmi))

  # Check if the number of levels is correct
  expect_equal(nlevels(discretized_data$age), 2)
  expect_equal(nlevels(discretized_data$bmi), 3)

  # Check that other columns are untouched
  expect_true(is.numeric(discretized_data$strokeha)) # It was numeric in test_data
})


# -----------------------------------------------------------------------------
# Smoke Test for the Main Generator Function
# -----------------------------------------------------------------------------

test_that("generate_synthetic_data runs without errors and returns a dataframe", {

  # This is a smoke test. It checks if the function runs to completion
  # without errors and returns an object of the expected class and dimensions.
  # It does not check the statistical properties of the output.

  # The function can be slow, so we use a very small dataset.
  # We need to create a dataframe with all the columns expected by the latent_connections list

  num_rows <- 100 # A bit larger to avoid issues with kmeans
  test_data_full <- data.frame(
    age = runif(num_rows, 40, 80),
    bmi = runif(num_rows, 20, 40),
    choleratio = runif(num_rows, 3, 7),
    sbp = runif(num_rows, 110, 160),
    sbps = runif(num_rows, 10, 30),
    strokeha = factor(sample(0:1, num_rows, replace = TRUE)),
    af = factor(sample(0:1, num_rows, replace = TRUE)),
    atyantip = factor(sample(0:1, num_rows, replace = TRUE)),
    steroid = factor(sample(0:1, num_rows, replace = TRUE)),
    impot = factor(sample(0:1, num_rows, replace = TRUE)),
    migr = factor(sample(0:1, num_rows, replace = TRUE)),
    ra = factor(sample(0:1, num_rows, replace = TRUE)),
    ckidney = factor(sample(0:1, num_rows, replace = TRUE)),
    semi = factor(sample(0:1, num_rows, replace = TRUE)),
    sle = factor(sample(0:1, num_rows, replace = TRUE)),
    treathyp = factor(sample(0:1, num_rows, replace = TRUE)),
    type1 = factor(sample(0:1, num_rows, replace = TRUE)),
    type2 = factor(sample(0:1, num_rows, replace = TRUE)),
    ethr = factor(sample(1:5, num_rows, replace = TRUE)),
    smoking = factor(sample(0:2, num_rows, replace = TRUE)),
    fh_cad = factor(sample(0:1, num_rows, replace = TRUE)),
    gender = factor(sample(c("M", "F"), num_rows, replace = TRUE)),
    region = factor(sample(1:10, num_rows, replace = TRUE))
  )

  var_names_before <- names(test_data_full)

  # Run the function and expect no errors
  synthetic_data <- NULL
  expect_no_error(synthetic_data <- generate_synthetic_data(test_data_full))

  # Check the output
  expect_true(is.data.frame(synthetic_data))

  # The number of rows can change due to the biological check (impotent females)
  # So we check if it's less than or equal to the original number of rows
  expect_lte(nrow(synthetic_data), nrow(test_data_full))

  # Check if latent variable columns are added
  num_latent_vars <- 6 # from config, but hardcoded here for test stability
  expect_equal(ncol(synthetic_data), ncol(test_data_full) + num_latent_vars)

  # Check that original columns are preserved (before latent vars are removed)
  expect_true(all(var_names_before %in% names(synthetic_data)))
})
