## Default 95% level

est = lm(yield ~ N + P + K, data = npk)

## Default 95% level
tidy_ri = tidy(ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L))
## Change level to 90 %
tidy_ri_90 = tidy(ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L, level = .90))

## Results are platform (actually: ptype) dependent
if (.Platform$OS.type != "windows") {
  tidy_ri_known = data.frame(
    term = 'N1',
    estimate = 5.616667,
    std.error = 0.007461833,
    p.value = 0.021,
    conf.low = 0.008726377,
    conf.high = 0.03327362
    )
  tidy_ri_90_known = data.frame(
    term = 'N1',
    estimate = 5.616667,
    std.error = 0.005813723,
    p.value = 0.021,
    conf.low = 0.01354941,
    conf.high = 0.02845059
  )
}

if (.Platform$OS.type == "windows") {
  tidy_ri_known = data.frame(
    term = 'N1',
    estimate = 5.616667,
    std.error = 0.008877507,
    p.value = 0.03,
    conf.low = 0.0153978,
    conf.high = 0.0446022
  )
  tidy_ri_90_known = data.frame(
    term = 'N1',
    estimate = 5.616667,
    std.error = 0.005813723,
    p.value = 0.021,
    conf.low = 0.01354941,
    conf.high = 0.02845059
  )
}

expect_equal(tidy_ri, tidy_ri_known, tolerance = 1e-6)

expect_equal(tidy_ri_90, tidy_ri_90_known, tolerance = 1e-6)

