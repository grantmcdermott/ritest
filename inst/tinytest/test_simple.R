#
# Simple RI. No strata (blocking) or clusters
#

# Use the example given in the `ri2` intro vignette, which is itself taken from
# Gerber and Green (2012, chpt 3)
table_2_2 = data.frame(Z = c(1, 0, 0, 0, 0, 0, 1),
                       Y = c(15, 15, 20, 20, 10, 15, 30))


est = lm(Y ~ Z, table_2_2)

#
# Forked processes
#

ritest_est = ritest(est, 'Z', reps = 1e3, seed = 42L, pcores = 2L)
ritest_est_fml = ritest(est, ~Z, reps = 1e3, seed = 42L, pcores = 2L)

# Test that character and formula results are the same
expect_equal(ritest_est, ritest_est_fml)

# Test p-value
expect_equal(as.numeric(ritest_est$pval), 0.373, tolerance = 1e-03)
# Test CI
expect_equal(as.numeric(ritest_est$ci), c(0.3316, 0.4144), tolerance = 1e-04)


#
# Sequential and PSOCK processes (will yield different results)
#

ritest_est_pcores1 = ritest(est, 'Z', reps = 1e3, seed = 42L, pcores = 1L)
ritest_est_seq = ritest(est, 'Z', reps = 1e3, seed = 42L, parallel = FALSE)

ritest_est_psock = ritest(est, 'Z', reps = 1e3, seed = 42L, pcores = 2L, ptype = 'psock')

# Test that PSOCK and sequential results are the same
expect_equal(ritest_est_pcores1, ritest_est_seq)

# Test p-value sequential
expect_equal(as.numeric(ritest_est_seq$pval), 0.383, tolerance = 1e-03)
# Test CI sequential
expect_equal(as.numeric(ritest_est_seq$ci), c(0.3414, 0.4246), tolerance = 1e-04)

# Test p-value psock
expect_equal(as.numeric(ritest_est_psock$pval), 0.384, tolerance = 1e-03)
# Test CI psock
expect_equal(as.numeric(ritest_est_psock$ci), c(0.3424, 0.4256), tolerance = 1e-04)
