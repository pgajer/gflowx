test_that("archive adapters reject objects outside their ownership boundary", {
  expect_error(
    lcor.with.rdgraph.posterior(list(), matrix(1, 1, 1)),
    "knn.riem.fit"
  )
  expect_error(
    lslope.test(list(), y = 1:3, z = 1:3, n.BB = 50, verbose = FALSE),
    "knn.riem.fit|requires optional package"
  )
})

test_that("subject neighborhood diagnostics are exported by gflowx", {
  expect_true("subject.neighborhood.stats" %in% getNamespaceExports("gflowx"))
})
