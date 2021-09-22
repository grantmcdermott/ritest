#
# Simple RI. No strata (blocking) or clusters
#

# Use the example given in the `ri2` intro vignette, which is itself taken from
# Gerber and Green (2012, chpt 3)
table_2_2 = data.frame(Z = c(1, 0, 0, 0, 0, 0, 1),
                       Y = c(15, 15, 20, 20, 10, 15, 30))


est = lm(Y ~ Z, table_2_2)

ritest_est = ritest(est, 'Z', reps = 1e3, seed = 42L, pcores = 2L)
ritest_est_fml = ritest(est, ~Z, reps = 1e3, seed = 42L, pcores = 2L)

# Test that character and formula results are the same
expect_equal(ritest_est, ritest_est_fml)

# Test p-value
expect_equal(as.numeric(ritest_est$pval), 0.373, tolerance = 1e-03)

# Test CI
expect_equal(as.numeric(ritest_est$ci), c(0.3316, 0.4144), tolerance = 1e-04)
