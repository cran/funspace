test_that("panel selection works", { # Testing find_nrow_ncol function
  expect_equal(find_nrow_ncol(group.vec = letters[1:6]), c(2, 3))
  expect_equal(find_nrow_ncol(group.vec = letters[1:4]), c(2, 2))
  expect_equal(find_nrow_ncol(group.vec = letters[1]), c(1, 1))
  expect_equal(find_nrow_ncol(group.vec = letters[1:12]), c(3, 4))
})

test_that("grid creation works", { # Testing grid_creation function
  expect_equal(nrow(grid_creation(list(1:100, 1:100), n_divisions = 100)) , 100 * 100)
  expect_equal(nrow(grid_creation(list(1:100, 1:100), n_divisions = 1)) , 1 * 1)
  expect_error(grid_creation(list(1:1, 1:1), n_divisions = 100))
})

