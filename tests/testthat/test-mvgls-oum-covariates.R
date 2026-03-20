make_oum_fixture <- function() {
  tree <- ape::rtree(8)
  tree$edge.length <- rep(1, nrow(tree$edge))
  states <- stats::setNames(rep(c("a", "b"), each = 4), tree$tip.label)
  simmap <- phytools::make.simmap(tree, x = states, model = "ER", nsim = 1, message = FALSE)

  Y <- cbind(
    trait1 = seq_len(8) / 4,
    trait2 = c(0.0, 0.5, 1.0, 1.5, -0.5, 0.0, 0.5, 1.0)
  )
  rownames(Y) <- tree$tip.label

  data <- data.frame(
    grp = factor(rep(c("g1", "g2"), each = 4)),
    x = seq(-1.4, 1.4, length.out = 8),
    row.names = tree$tip.label
  )

  list(tree = simmap, Y = Y, data = data)
}

test_that("OUM retains covariates for downstream regression helpers", {
  set.seed(1)
  fixture <- make_oum_fixture()

  fit <- suppressWarnings(
    mvgls(
      fixture$Y ~ grp + x,
      data = fixture$data,
      tree = fixture$tree,
      model = "OUM",
      method = "LL"
    )
  )

  coef_names <- rownames(coef(fit))
  expect_true("x" %in% coef_names)
  expect_true(any(grepl("^grp", coef_names)))
  expect_equal(nrow(coef(fit)), ncol(fit$variables$regimes) + 2L)

  predicted <- predict(fit, newdata = fixture$data)
  expect_equal(dim(predicted), dim(fixture$Y))

  node_data <- data.frame(
    grp = factor(rep("g1", ape::Nnode(fixture$tree)), levels = levels(fixture$data$grp)),
    x = seq(-0.5, 0.5, length.out = ape::Nnode(fixture$tree))
  )
  ancestral_fit <- ancestral(fit, newdata = node_data)
  expect_equal(dim(ancestral_fit), c(ape::Nnode(fixture$tree), ncol(fixture$Y)))

  mv_test <- manova.gls(fit, type = "II", test = "Pillai")
  expect_s3_class(mv_test, "manova.mvgls")
  expect_true(all(c("grp", "x") %in% mv_test$terms))
})

test_that("OUM information criteria work for stored formula and data calls", {
  set.seed(2)
  fixture <- make_oum_fixture()
  Y <- fixture$Y
  predictors <- fixture$data
  model_formula <- Y ~ grp + x

  fit <- suppressWarnings(
    mvgls(
      model_formula,
      data = predictors,
      tree = fixture$tree,
      model = "OUM",
      method = "LL"
    )
  )

  gic <- GIC(fit)
  expect_s3_class(gic, "gic.mvgls")
  expect_true(is.finite(gic$GIC))

  eic <- suppressWarnings(EIC(fit, nboot = 1L, nbcores = 1L))
  expect_s3_class(eic, "eic.mvgls")
  expect_true(is.finite(eic$EIC))
  expect_true(is.na(eic$se) || is.finite(eic$se))
})

test_that("OUM print output separates regime optima from regression effects", {
  set.seed(3)
  fixture <- make_oum_fixture()

  fit <- suppressWarnings(
    mvgls(
      fixture$Y ~ grp + x,
      data = fixture$data,
      tree = fixture$tree,
      model = "OUM",
      method = "LL"
    )
  )

  printed <- capture.output(print(fit))
  expect_true(any(grepl("Regime optima \\(theta\\):", printed)))
  expect_true(any(grepl("Regression effects \\(beta\\):", printed)))

  summary_printed <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Regime optima \\(theta\\):", summary_printed)))
  expect_true(any(grepl("Regression effects \\(beta\\):", summary_printed)))
})

test_that("OUM works on the phyllostomid example dataset with a covariate", {
  set.seed(123)
  data("phyllostomid", package = "mvMORPH")

  states <- phyllostomid$grp2
  names(states) <- names(phyllostomid$grp2)
  tree_sim <- suppressMessages(
    phytools::make.simmap(phyllostomid$tree, x = states, model = "ER", nsim = 1, message = FALSE)
  )

  Y <- phyllostomid$mandible[, 1:6, drop = FALSE]
  size <- as.numeric(scale(phyllostomid$mandible[, 73]))
  dat <- list(Y = Y, size = size)

  fit <- suppressWarnings(
    mvgls(
      Y ~ size,
      data = dat,
      tree = tree_sim,
      model = "OUM",
      method = "LL",
      REML = FALSE
    )
  )

  expect_true("size" %in% rownames(coef(fit)))
  expect_equal(dim(fit$variables$X), c(nrow(Y), 4L))

  predicted <- predict(fit, newdata = data.frame(size = size, row.names = tree_sim$tip.label))
  expect_equal(dim(predicted), dim(Y))

  mv_test <- manova.gls(fit, type = "II", test = "Pillai")
  expect_identical(mv_test$terms, "size")

  gic <- GIC(fit)
  eic <- suppressWarnings(EIC(fit, nboot = 1L, nbcores = 1L))
  expect_true(is.finite(gic$GIC))
  expect_true(is.finite(eic$EIC))
})
