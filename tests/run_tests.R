# This script runs all the unit tests for the project.

# Load the testthat library
library(testthat)

# Run all tests in the 'testthat' directory
test_dir("tests/testthat/", reporter = "summary")
