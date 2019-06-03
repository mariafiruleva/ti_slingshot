test_that("ti_slingshot works", {
  dataset <- source(system.file("example.sh", package = "tislingshot"))$value

  model <- dynwrap::infer_trajectory(dataset, tislingshot::ti_slingshot(), verbose = T)

  expect_is(model, "dynwrap::with_trajectory")
})
