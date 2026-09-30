# xgboost backend: fit an xgboost model to the lagged features and predict the next value
backend_xgboost <- function() {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("install.packages('xgboost')")
    backend_lagmodel("xgboost", function(tab) {
        x <- as.matrix(tab[, setdiff(names(tab), "y")])
        y <- tab$y
        fit <- xgboost::xgboost(
            data = x, label = y, nrounds = 150, max_depth = 4, eta = 0.1,
            objective = "reg:quantileerror", quantile_alpha = LEVELS, verbose = 0
        )
        function(newdata) matrix(predict(fit, as.matrix(newdata), reshape = TRUE), ncol = length(LEVELS))
    })
}