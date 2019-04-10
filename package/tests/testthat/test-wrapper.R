test_that("ti_slingshot works", {
  dataset <- dynutils::read_h5(system.file("example.h5", package = "tislingshot"))

  model <- dynwrap::infer_trajectory(dataset, tislingshot::ti_slingshot(), verbose = T)

  expect_is(model, "dynwrap::with_trajectory")
})
