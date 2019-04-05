#!/usr/local/bin/Rscript

requireNamespace("dyncli", quietly = TRUE)
task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)

requireNamespace("princurve", quietly = TRUE)
requireNamespace("cluster", quietly = TRUE)
requireNamespace("irlba", quietly = TRUE)

suppressWarnings(library(slingshot, warn.conflicts = FALSE))

#####################################
###           LOAD DATA           ###
#####################################

#' @example
#' system("Rscript ti_slingshot/example.sh /tmp/test.h5")
#' task <- dyncli::main(
#'   args = "--dataset /tmp/test.h5 --output /tmp/output.h5" %>% strsplit(" ") %>% first(),
#'   definition_location = "ti_slingshot/definition.yml"
#' )

parameters <- task$parameters
expression <- task$expression
start_id <- task$priors$start_id
end_id <- task$priors$end_id

#####################################
###        INFER TRAJECTORY       ###
#####################################

#   ____________________________________________________________________________
#   Preprocessing                                                           ####

start_cell <- if (!is.null(start_id)) { sample(start_id, 1) } else { NULL }

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

#   ____________________________________________________________________________
#   Dimensionality reduction                                                ####
pca <- irlba::prcomp_irlba(expression, n = 20)

# this code is adapted from the expermclust() function in TSCAN
# the only difference is in how PCA is performed
# (they specify scale. = TRUE and we leave it as FALSE)
x <- 1:20
optpoint1 <- which.min(sapply(2:10, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(pca$sdev[1:20] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
}))

# this is a simple method for finding the "elbow" of a curve, from
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
x <- cbind(1:20, pca$sdev[1:20])
line <- x[c(1, nrow(x)),]
proj <- princurve::project_to_curve(x, line)
optpoint2 <- which.max(proj$dist_ind)-1

# we will take more than 3 PCs only if both methods recommend it
optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
rd <- pca$x[, seq_len(optpoint)]
rownames(rd) <- rownames(expression)

#   ____________________________________________________________________________
#   Clustering                                                              ####
# max clusters equal to number of cells
max_clusters <- min(nrow(expression)-1, 10)

clusterings <- lapply(3:max_clusters, function(K){
  cluster::pam(rd, K) # we generally prefer PAM as a more robust alternative to k-means
})

# take one more than the optimal number of clusters based on average silhouette width
# (max of 10; the extra cluster improves flexibility when learning the topology,
# silhouette width tends to pick too few clusters, otherwise)
wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
labels <- clusterings[[min(c(wh.cl, 8))]]$clustering

start.clus <-
  if(!is.null(start_cell)) {
    labels[[start_cell]]
  } else {
    NULL
  }
end.clus <-
  if(!is.null(end_id)) {
    unique(labels[end_id])
  } else {
    NULL
  }

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
sds <- slingshot(
  rd,
  labels,
  start.clus = start.clus,
  end.clus = end.clus,
  shrink = parameters$shrink,
  reweight = parameters$reweight,
  reassign = parameters$reassign,
  thresh = parameters$thresh,
  maxit = parameters$maxit,
  stretch = parameters$stretch,
  smoother = parameters$smoother,
  shrink.method = parameters$shrink.method
)

start_cell <- apply(slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
start.clus <- labels[[start_cell]]

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Create output                                                           ####

# collect milestone network
lineages <- slingLineages(sds)
lineage_ctrl <- slingParams(sds)

cluster_network <- lineages %>%
  map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
  unique() %>%
  mutate(
    length = lineage_ctrl$dist[cbind(from, to)],
    directed = TRUE
  )

# collect dimred
dimred <- reducedDim(sds)

# collect clusters
cluster <- clusterLabels(sds)

# collect progressions
adj <- slingAdjacency(sds)
lin_assign <- apply(slingCurveWeights(sds), 1, which.max)

progressions <- map_df(seq_along(lineages), function(l) {
  ind <- lin_assign == l
  lin <- lineages[[l]]
  pst.full <- slingPseudotime(sds, na = FALSE)[,l]
  pst <- pst.full[ind]
  means <- sapply(lin, function(clID){
    stats::weighted.mean(pst.full, cluster[,clID])
  })
  non_ends <- means[-c(1,length(means))]
  edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
  from.l <- lineages[[l]][edgeID.l]
  to.l <- lineages[[l]][edgeID.l + 1]
  m.from <- means[from.l]
  m.to <- means[to.l]

  pct <- (pst - m.from) / (m.to - m.from)
  pct[pct < 0] <- 0
  pct[pct > 1] <- 1

  tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
})

#   ____________________________________________________________________________
#   Save output                                                             ####
output <-
  wrap_data(
    cell_ids = rownames(expression)
  ) %>%
  add_trajectory(
    milestone_network = cluster_network,
    progressions = progressions
  ) %>%
  add_dimred(
    dimred = dimred
  ) %>%
  add_timings(checkpoints)

dyncli::write_output(output, task$output)

