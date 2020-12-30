#!/usr/bin/env Rscript

requireNamespace("dyncli", quietly = TRUE)
task <- dyncli::main()

library(tislingshot, warn.conflicts = FALSE)

output <- tislingshot::run_fun(
  expression = task$expression,
  priors = task$priors,
  parameters = task$parameters,
  seed = task$seed,
  verbose = task$verbose
)

output$sling_out <- list()
output$sling_out$pseudotime <- slingshot::slingPseudotime(sds) %>% as.data.frame()
output$sling_out$cell_weights <- slingshot::slingCurveWeights(sds) %>% as.data.frame()

dyncli::write_output(output, task$output)
