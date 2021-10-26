#
# RI with strata and clusters
#

# Use the bundled 'colombia' dataset
data("colombia")

co_est =
  fixest::feols(
    dayscorab ~ b_treat + b_dayscorab + miss_b_dayscorab + round2 + round3 | b_pair,
    vcov = ~b_block, data = colombia
    )

## Default parallel (with 2 cores)
co_ri = ritest(co_est, ~b_treat, cluster=~b_block, strata=~b_pair,
               reps = 1e3, seed = 546L, pcores = 2L)
## Sequential
co_ri_seq = ritest(co_est, ~b_treat, cluster=~b_block, strata=~b_pair,
                   reps = 1e3, seed = 546L, parallel = FALSE)
## Parallel, but no stacking (will be different)
co_ri_nostack = ritest(co_est, ~b_treat, cluster=~b_block, strata=~b_pair,
                       reps = 1e3, seed = 546L, pcores = 2L, stack = FALSE)
## Sequential and no stacking
co_ri_seq_nostack = ritest(co_est, ~b_treat, cluster=~b_block, strata=~b_pair,
                           reps = 1e3, seed = 546L, parallel = FALSE, stack = FALSE)

# Test that default and sequential results are the same
expect_equal(co_ri, co_ri_seq)
# Test that default and sequential+non-stacked results are the same
expect_equal(co_ri, co_ri_seq_nostack)

# Test count
expect_equal(co_ri$count, 102)
# Test SE
expect_equal(co_ri$se, 0.01575008)

# Test count (no stacking)
expect_equal(co_ri_nostack$count, 119)
# Test SE (no stacking)
expect_equal(co_ri_nostack$se, 0.01685023, tolerance = 1e-6)
