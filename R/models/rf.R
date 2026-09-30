LAGS <- c(1, 2, 3, 5, 7, 10, 14, 21, 28, 35, 49, 63, 90)
ROLL <- c(7, 14, 30)

# build a lag table for a given time series, with specified lags and rolling windows
build_lag_table <- function(x, lags = LAGS, roll = ROLL) {
    n <- length(x)
    maxlag <- max(lags, roll)
    idx <- (maxlag + 1):(n - 1)
    if (length(idx) < 60) {
        return(NULL)
    }
    feats <- lapply(lags, function(L) x[idx - L])
    names(feats) <- paste0("lag", lags)
    for (w in roll) {
        feats[[paste0("ma", w)]] <- vapply(idx, function(t) mean(x[(t - w + 1):t]), numeric(1))
        feats[[paste0("sd", w)]] <- vapply(idx, function(t) sd(x[(t - w + 1):t]), numeric(1))
    }
    df <- as.data.frame(feats)
    df$y <- x[idx + 1]
    df
}
# compute the lagged features for a single row (the last row of a time series)
lag_state <- function(x, lags = LAGS, roll = ROLL) {
    n <- length(x)
    row <- as.list(setNames(x[n - lags + 1], paste0("lag", lags)))
    for (w in roll) {
        row[[paste0("ma", w)]] <- mean(x[(n - w + 1):n])
        row[[paste0("sd", w)]] <- sd(x[(n - w + 1):n])
    }
    as.data.frame(row)
}
# sample from a matrix of quantiles, using linear interpolation for the inner quantiles and normal approximation for the tails
sample_row_quantiles <- function(Q, u) {
    n <- nrow(Q)
    out <- numeric(n)
    inner <- u >= 0.1 & u <= 0.9
    lo <- u < 0.1
    hi <- u > 0.9
    s_lo <- pmax((Q[, 5] - Q[, 1]) / Z90, 1e-6)
    s_hi <- pmax((Q[, 9] - Q[, 5]) / Z90, 1e-6)
    if (any(inner)) {
        out[inner] <- vapply(
            which(inner),
            function(i) approx(LEVELS, Q[i, ], xout = u[i], rule = 2)$y, numeric(1)
        )
    }
    out[lo] <- Q[lo, 1] + s_lo[lo] * (qnorm(u[lo]) - qnorm(0.1))
    out[hi] <- Q[hi, 9] + s_hi[hi] * (qnorm(u[hi]) - qnorm(0.9))
    out
}
# simulate multiple paths from a fitted model, using the lagged features and the step function
simulate_paths <- function(x0, step_fn, horizon, n) {
    mat <- matrix(rep(x0, n), nrow = n, byrow = TRUE)
    out <- matrix(NA_real_, n, horizon)
    u_grid <- (seq_len(n) - 0.5) / n
    for (h in seq_len(horizon)) {
        state <- do.call(rbind, lapply(seq_len(n), function(i) lag_state(mat[i, ])))
        Q <- step_fn(state)
        out[, h] <- sample_row_quantiles(Q, sample(u_grid))
        mat <- cbind(mat, out[, h])
    }
    out
}
# fit a lagged model and return a function that predicts the next value given the lagged features
backend_lagmodel <- function(name, fit_and_predict) {
    list(name = name, predict = function(ctxs, horizon) {
        out <- array(NA_real_, c(length(ctxs), horizon, 9))
        for (i in seq_along(ctxs)) {
            x <- ctxs[[i]]
            tab <- build_lag_table(x)
            if (is.null(tab)) {
                out[i, , ] <- x[length(x)]
                next
            }
            step_fn <- fit_and_predict(tab)
            paths <- simulate_paths(x, step_fn, horizon, 199)
            out[i, , ] <- t(apply(paths, 2, quantile, probs = LEVELS, names = FALSE))
        }
        out
    })
}
# random forest backend: fit a random forest to the lagged features and predict the next value
backend_rf <- function() {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("install.packages('ranger')")
    backend_lagmodel("rf", function(tab) {
        fit <- ranger::ranger(y ~ ., data = tab, quantreg = TRUE, num.trees = 200)
        function(newdata) predict(fit, data = newdata, type = "quantiles", quantiles = LEVELS)$predictions
    })
}