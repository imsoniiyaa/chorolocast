DIRECT_H <- c(1, 3, 7, 10, 14, 18, 21, 25, 28, 32, 35)

direct_features <- function(x, ts, lags = LAGS, roll = ROLL) {
    feats <- lapply(lags, function(L) x[ts - L + 1])
    names(feats) <- paste0("lag", lags)
    for (w in roll) {
        feats[[paste0("ma", w)]] <- vapply(ts, function(t) mean(x[(t - w + 1):t]), numeric(1))
        feats[[paste0("sd", w)]] <- vapply(ts, function(t) sd(x[(t - w + 1):t]), numeric(1))
    }
    as.data.frame(feats)
}

backend_rf_direct <- function(num_trees = 300) {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("install.packages('ranger')")
    seed <- as.integer(Sys.getenv("SEED", "42"))
    K <- length(LEVELS)
    maxlag <- max(c(LAGS, ROLL))

    list(name = "rf_direct", predict = function(ctxs, horizon) {
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
                tab$y <- x[ts + h] - x[ts]

                fit <- ranger::ranger(
                    y ~ .,
                    data = tab, quantreg = TRUE,
                    num.trees = num_trees, seed = seed
                )
                Qg[g, ] <- x[n] + predict(fit, data = newrow, type = "quantiles", quantiles = LEVELS)$predictions[1, ]
            }
            for (k in seq_len(K)) {
                out[i, , k] <- approx(grid, Qg[, k], xout = seq_len(horizon), rule = 2)$y
            }
        }
        out
    })
}
