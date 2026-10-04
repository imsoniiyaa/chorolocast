backend_arima <- function(period = 365.25, Ks = c(0, 2, 4), min_cycles = 1.2,
                          use_do = TRUE,
                          gamma = as.numeric(Sys.getenv("DO_GAMMA", "0.5")),
                          hgrid = c(1, 3, 7, 14, 21, 28, 35)) {
    if (!requireNamespace("forecast", quietly = TRUE)) stop("install.packages('forecast')")

    fourier_terms <- function(t, K) {
        if (K == 0) {
            return(NULL)
        }
        m <- do.call(cbind, lapply(1:K, function(k) {
            cbind(sin(2 * pi * k * t / period), cos(2 * pi * k * t / period))
        }))
        colnames(m) <- paste0(c("S", "C"), rep(1:K, each = 2))
        m
    }

    fit_K <- function(z, K) {
        xr <- fourier_terms(seq_along(z), K)
        fit <- tryCatch(
            forecast::auto.arima(z, xreg = xr, seasonal = FALSE, stepwise = TRUE),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            return(NULL)
        }
        fit$call$xreg <- xr
        list(fit = fit, K = K, aicc = fit$aicc)
    }

    # context 안에서 h별 DOsat 계수 추정
    do_beta <- function(y, dosat, hs) {
        n <- length(y)
        if (is.null(dosat) || sum(!is.na(dosat)) < 60) {
            return(NULL)
        }
        d <- approx(seq_len(n), dosat, xout = seq_len(n), rule = 2)$y
        if (!is.finite(sd(d)) || sd(d) < 1e-8) {
            return(NULL)
        }
        dz <- (d - mean(d)) / sd(d)
        ctrl <- as.numeric(y - stats::filter(y, rep(1 / 30, 30), sides = 1))
        beta <- numeric(length(hs))
        for (g in seq_along(hs)) {
            h <- hs[g]
            ts <- 30:(n - h)
            if (length(ts) < 60) next
            tab <- data.frame(dy = y[ts + h] - y[ts], dz = dz[ts], ctrl = ctrl[ts])
            fit <- tryCatch(lm(dy ~ dz + ctrl, data = tab), error = function(e) NULL)
            if (!is.null(fit) && is.finite(coef(fit)["dz"])) beta[g] <- unname(coef(fit)["dz"])
        }
        list(beta = beta, dz_last = dz[n])
    }

    list(name = "arima", predict = function(ctxs, horizon) {
        out <- array(NA_real_, c(length(ctxs), horizon, length(LEVELS)))

        for (i in seq_along(ctxs)) {
            x <- if (is.list(ctxs[[i]])) ctxs[[i]]$y else ctxs[[i]]
            n <- length(x)

            Kcand <- if (n >= min_cycles * period) Ks else 0
            fits <- Filter(Negate(is.null), lapply(Kcand, function(K) fit_K(x, K)))

            if (length(fits) == 0) {
                fit <- forecast::auto.arima(x, seasonal = FALSE)
                K <- 0
            } else {
                best <- fits[[which.min(sapply(fits, `[[`, "aicc"))]]
                fit <- best$fit
                K <- best$K
            }

            newxreg <- fourier_terms(n + seq_len(horizon), K)
            fc <- forecast::forecast(fit, h = horizon, xreg = newxreg, level = 95)
            se <- pmax((fc$upper[, 1] - fc$mean) / qnorm(0.975), 1e-6)

            ramp <- pmin(1, pmax(0, (seq_len(horizon) - 3) / 11))
            corr <- rep(0, horizon)
            if (use_do && gamma > 0 && is.list(ctxs[[i]])) {
                hs <- sort(unique(c(hgrid[hgrid <= horizon], horizon)))
                b <- do_beta(x, ctxs[[i]]$dosat, hs)
                if (!is.null(b)) {
                    bh <- if (length(hs) < 2) {
                        rep(b$beta[1], horizon)
                    } else {
                        approx(hs, b$beta, xout = seq_len(horizon), rule = 2)$y
                    }
                    ramp <- pmin(1, pmax(0, (seq_len(horizon) - 3) / 11))
                    corr <- pmin(pmax(gamma * bh * b$dz_last, -0.5), 0.5) * ramp
                }
            }

            for (h in seq_len(horizon)) {
                out[i, h, ] <- as.numeric(fc$mean[h]) + corr[h] + se[h] * qnorm(LEVELS)
            }
        }
        out
    })
}
