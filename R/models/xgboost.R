# xgboost backend
backend_xgboost <- function() {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("install.packages('xgboost')")
    seed <- as.integer(Sys.getenv("SEED", "42"))
    K <- length(LEVELS)
    backend_lagmodel("xgboost", function(tab) {
        x <- as.matrix(tab[, setdiff(names(tab), "y")])
        dtrain <- xgboost::xgb.DMatrix(data = x, label = tab$y)
        fit <- xgboost::xgb.train(
            params = list(
                objective = "reg:quantileerror",
                quantile_alpha = LEVELS,
                max_depth = 4,
                eta = 0.03,
                subsample = 0.8,
                colsample_bytree = 0.8,
                seed = seed
            ),
            data = dtrain,
            nrounds = 300,
            verbose = 0
        )
        function(newdata) {
            p <- predict(fit, xgboost::xgb.DMatrix(as.matrix(newdata[, colnames(x), drop = FALSE])))
            if (is.null(dim(p))) p <- matrix(p, ncol = K)
            t(apply(p, 1, sort)) # quantile crossing 방지
        }
    })
}