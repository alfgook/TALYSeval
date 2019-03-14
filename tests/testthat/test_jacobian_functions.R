context("Jacobian functions")
library(data.table)

needsDt <- data.table(IDX = 1:4,
                      PROJECTILE = c("N", "N", "N", "N"),
                      ELEMENT = c("FE", "FE", "FE", "FE"),
                      MASS = c(54L, 54L, 56L, 56L),
                      REAC = rep("CS/TOT", 4),
                      L1 = c(1, 2, 3, 4), L2 = 0, L3 = 0)

paramDt <- data.table(IDX = 1:6,
                      PROJECTILE = rep("N", 6),
                      ELEMENT = rep("FE", 6),
                      MASS = c(rep(54L, 3), rep(56L, 3)),
                      PARNAME = rep(c("a", "b", "energy"), 2),
                      PARVAL = list(1, 2, seq(1,10,by=1),
                                    3, 4, seq(1,10,by=1)),
                      ADJUSTABLE = c(TRUE, TRUE, FALSE))


talysMockup <- function(inputs, outspecs) {

  mapply(function(input, outspec) {

    outspec <- copy(outspec)
    outspec[, V1 := L1 * input$a^2 + L1 * input$b]
  }, input = inputs, outspec = outspecs, SIMPLIFY = FALSE)
}


test_that("sensitivity coefficients are determined correctly", {

  variationDt <- createInputsForJacobian(paramDt, needsDt, eps = 1e-6)
  outspecs <- talysMockup(variationDt$inputs, variationDt$outspecs)
  variationDt$outspecs <- outspecs
  jacRes <- computeJacobian(variationDt)
  # prepare results
  setkey(jacRes, SRCIDX, DSTIDX)
  setkey(paramDt, IDX)
  setkey(needsDt, IDX)
  tmpDt <- cbind(paramDt[J(jacRes$SRCIDX),], needsDt[J(jacRes$DSTIDX),], estS = jacRes$S)

  tmpDt[, trueS := {
    if (.BY[["PARNAME"]] == "a")
      2 * unlist(PARVAL) * L1
    else if (.BY[["PARNAME"]] == "b")
      L1
  }, by="PARNAME"]

  expect_equal(tmpDt$estS, tmpDt$trueS, tolerance=1e-4)
})




