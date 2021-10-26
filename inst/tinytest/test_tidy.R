## Default 95% level

est = lm(yield ~ N + P + K, data = npk)

tidy_ri = tidy(ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L))

tidy_ri_known = data.frame(
  term = 'N1',
  estimate = 5.616667,
  std.error = 0.007461833,
  p.value = 0.021,
  conf.low = 0.008726377,
  conf.high = 0.03327362
  )

expect_equal(tidy_ri, tidy_ri_known, tolerance = 1e-6)

## Change level to 90 %

tidy_ri_90 = tidy(ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L, level = .90))

tidy_ri_90_known = data.frame(
  term = 'N1',
  estimate = 5.616667,
  std.error = 0.005813723,
  p.value = 0.021,
  conf.low = 0.01354941,
  conf.high = 0.02845059
)

expect_equal(tidy_ri_90, tidy_ri_90_known, tolerance = 1e-6)

