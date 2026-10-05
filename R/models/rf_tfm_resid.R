b <- backend_timesfm(512)
get_timesfm_quantiles <- function(x, horizon) {
    arr <- b$predict(
        list(x),
        horizon
    )

    as.numeric(arr[1, 1, ])
}


direct_features <- function(x, ts, lags = LAGS, roll = ROLL) {
    feats <- lapply(lags, function(L) x[ts - L + 1])
    names(feats) <- paste0("lag", lags)
    for (w in roll) {
        feats[[paste0("ma", w)]] <- vapply(ts, function(t) mean(x[(t - w + 1):t]), numeric(1))
        feats[[paste0("sd", w)]] <- vapply(ts, function(t) sd(x[(t - w + 1):t]), numeric(1))
    }
    as.data.frame(feats)
}
backend_rf_tfm_resid <- function(num_trees = 300) {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("install.packages('ranger')")
    seed <- as.integer(Sys.getenv("SEED", "42"))
    K <- length(LEVELS)
    maxlag <- max(c(LAGS, ROLL))

    list(name = "rf_tfm_resid", predict = function(ctxs, horizon) {
        out <- array(NA_real_, c(length(ctxs), horizon, K))
        grid <- sort(unique(c(DIRECT_H[DIRECT_H <= horizon], horizon)))

        for (i in seq_along(ctxs)) {
            x <- if (is.list(ctxs[[i]])) ctxs[[i]]$y else ctxs[[i]]
            n <- length(x)
            if (n - max(grid) - maxlag + 1 < 60) {
                out[i, , ] <- x[n]
                next
            }
            ts_all <- maxlag:n
            F_all <- direct_features(x, ts_all)
            newrow <- F_all[nrow(F_all), , drop = FALSE]

            Qg <- matrix(NA_real_, length(grid), K)
            for (g in seq_along(grid)) {
                h <- grid[g]
                ts <- maxlag:(n - h)
                tab <- F_all[seq_along(ts), , drop = FALSE]
                tfm_train <- numeric(length(ts))
                for (j in seq_along(ts)) {
                    ctx_train <- x[1:ts[j]]

                    Q <- get_timesfm_quantiles(
                        ctx_train,
                        h
                    )

                    tfm_train[j] <- Q[5]
                }
                tab$y <- x[ts + h] - tfm_train

                fit <- ranger::ranger(
                    y ~ .,
                    data = tab,
                    quantreg = TRUE,
                    num.trees = num_trees,
                    seed = seed
                )

                tfm_future <- get_timesfm_quantiles(
                    x,
                    h
                )

                resid_q <- predict(
                    fit,
                    data = newrow,
                    type = "quantiles",
                    quantiles = LEVELS
                )$predictions[1, ]

                Qg[g, ] <- tfm_future + resid_q
            }
            for (k in seq_len(K)) {
                out[i, , k] <- approx(grid, Qg[, k], xout = seq_len(horizon), rule = 2)$y
            }
        }
        out
    })
}
