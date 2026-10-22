# Exercise the fitting guard without running FUSION's top-level pipeline.
expressions = parse("pheast/scripts/fusion_twas/FUSION.compute_weights.R")
for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("=")) &&
        identical(expr[[2]], as.name("weights.enet"))) eval(expr)
}

cleanup = function() cleaned <<- TRUE
q = function(...) stop(structure(list(message="scientific skip"),
                                  class=c("scientific_skip", "error", "condition")))
expect.skip = function(x, diagnostic) {
    cleaned <<- FALSE
    output = capture.output(result <- tryCatch(
        weights.enet(x, matrix(seq_len(nrow(x)))),
        scientific_skip=function(e) e))
    stopifnot(inherits(result, "scientific_skip"), cleaned,
              any(grepl(diagnostic, output, fixed=TRUE)))
}
cv.glmnet = function(...) stop("glmnet must not be called for a degenerate split")
expect.skip(matrix(0, 41, 3), "full training set")
# One carrier in an outer training set becomes zero in one inner split.
x = matrix(rep(c(1, rep(0, 40)), 3), ncol=3)
set.seed(1)
expect.skip(x, "inner training fold")

# A second informative predictor keeps the split trainable. Check that the
# exact validated folds reach glmnet and filtered weights retain SNP order.
x = cbind(x, seq_len(41), 0)
cv.glmnet = function(x, y, alpha, nfolds, foldid, intercept, standardize) {
    stopifnot(is.matrix(x), ncol(x) == 4, nfolds == 5,
              length(foldid) == nrow(x))
    for (fold in unique(foldid)) {
        stopifnot(any(apply(x[foldid != fold,,drop=FALSE], 2, sd) > 0))
    }
    structure(list(), class="test_enet")
}
coef.test_enet = function(object, s, ...) matrix(0:4)
stopifnot(identical(as.vector(weights.enet(x, matrix(1:41))), c(1:4, 0)))

# Unexpected fitting errors propagate, rather than becoming scientific skips.
cv.glmnet = function(...) stop("unexpected fitting error")
result = tryCatch(weights.enet(x, matrix(1:41)), error=function(e) e)
stopifnot(identical(conditionMessage(result), "unexpected fitting error"))
cat("Elastic-net CV guard checks passed\n")
