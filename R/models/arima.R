# ARIMA backend: fit an ARIMA model to the time series and predict the next values,
# using Fourier terms for seasonality if enough data is available
backend_arima <- function() {
    if (!requireNamespace("forecast", quietly = TRUE)) stop("install.packages('forecast')")
    order_cache <- NULL
    fourier <- function(t, K = 4, period = 365.25) {
        m <- do.call(cbind, lapply(1:K, function(k) cbind(sin(2 * pi * k * t / period), cos(2 * pi * k * t / period))))
        colnames(m) <- paste0(c("S", "C"), rep(1:K, each = 2))
        m
    }
    list(name = "arima", predict = function(ctxs, horizon) {
        out <- array(NA_real_, c(length(ctxs), horizon, 9))
        for (i in seq_along(ctxs)) {
            x <- ctxs[[i]]
            n <- length(x)
            xreg <- if (n >= 2 * 365) fourier(seq_len(n)) else NULL
            fit <- tryCatch(
                {
                    if (is.null(order_cache)) {
                        f <- forecast::auto.arima(x, xreg = xreg, seasonal = FALSE, stepwise = TRUE, approximation = TRUE)
                        order_cache <<- as.integer(forecast::arimaorder(f)[1:3])
                        f
                    } else {
                        forecast::Arima(x, order = order_cache, xreg = xreg)
                    }
                },
                error = function(e) forecast::auto.arima(x, seasonal = FALSE)
            )
            newxreg <- if (!is.null(xreg)) fourier(n + seq_len(horizon)) else NULL
            fc <- forecast::forecast(fit, h = horizon, xreg = newxreg, level = 95)
            se <- pmax((fc$upper[, 1] - fc$mean) / qnorm(0.975), 1e-6)
            for (h in seq_len(horizon)) out[i, h, ] <- as.numeric(fc$mean[h]) + se[h] * qnorm(LEVELS)
        }
        out
    })
}