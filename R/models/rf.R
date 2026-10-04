LAGS <- c(1:7, 14, 30, 60, 90, 180, 365)
ROLL <- c(7, 14, 30)

if (!exists("LEVELS")) LEVELS <- seq(0.1, 0.9, by = 0.1)
if (!exists("Z90")) Z90 <- qnorm(0.9)

build_lag_table <- function(
  x,
  do = NULL,
  dosat = NULL,
  temp = NULL,
  airtemp = NULL,
  lags = LAGS,
  roll = ROLL
) {
    n <- length(x)
    maxlag <- max(c(lags, roll))
    idx <- (maxlag + 1):(n - 1)
    if (length(idx) < 60) {
        return(NULL)
    }
    feats <- lapply(lags, function(L) x[idx + 1 - L])
    names(feats) <- paste0("lag", lags)
    for (w in roll) {
        feats[[paste0("ma", w)]] <- vapply(idx, function(t) mean(x[(t - w + 1):t]), numeric(1))
        feats[[paste0("sd", w)]] <- vapply(idx, function(t) sd(x[(t - w + 1):t]), numeric(1))
    }
    if (!is.null(do)) {
        feats$DO <- do[idx]
    }

    if (!is.null(dosat)) {
        feats$DOsat <- dosat[idx]
    }

    if (!is.null(temp)) {
        feats$Temp <- temp[idx]
    }

    if (!is.null(airtemp)) {
        feats$AirTemp <- airtemp[idx]
    }
    df <- as.data.frame(feats)
    df$y <- x[idx + 1]
    df
}

# mat의 각 행(path)에서 마지막 시점 기준 feature를 한 번에 계산
lag_state_mat <- function(
  mat,
  do_last,
  dosat_last,
  temp_last,
  airtemp_last,
  lags = LAGS,
  roll = ROLL
) {
    n <- ncol(mat)
    out <- lapply(lags, function(L) mat[, n - L + 1])
    names(out) <- paste0("lag", lags)
    for (w in roll) {
        win <- mat[, (n - w + 1):n, drop = FALSE]
        out[[paste0("ma", w)]] <- rowMeans(win)
        out[[paste0("sd", w)]] <- apply(win, 1, sd)
    }
    out$DO <- rep(do_last, nrow(mat))
    out$DOsat <- rep(dosat_last, nrow(mat))
    out$Temp <- rep(temp_last, nrow(mat))
    out$AirTemp <- rep(airtemp_last, nrow(mat))

    as.data.frame(out)
}

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
            function(i) approx(LEVELS, Q[i, ], xout = u[i], rule = 2)$y,
            numeric(1)
        )
    }
    out[lo] <- Q[lo, 1] + s_lo[lo] * (qnorm(u[lo]) - qnorm(0.1))
    out[hi] <- Q[hi, 9] + s_hi[hi] * (qnorm(u[hi]) - qnorm(0.9))
    out
}

simulate_paths <- function(
  x0,
  step_fn,
  horizon,
  n,
  do_last,
  dosat_last,
  temp_last,
  airtemp_last
) {
    mat <- matrix(rep(x0, n), nrow = n, byrow = TRUE)
    out <- matrix(NA_real_, n, horizon)
    u_grid <- (seq_len(n) - 0.5) / n
    for (h in seq_len(horizon)) {
        state <- lag_state_mat(
            mat,
            do_last,
            dosat_last,
            temp_last,
            airtemp_last
        )
        Q <- step_fn(state)
        out[, h] <- sample_row_quantiles(Q, sample(u_grid))
        mat <- cbind(mat, out[, h])
    }
    out
}

backend_lagmodel <- function(name, fit_and_predict) {
    list(
        name = name,
        predict = function(ctxs, horizon) {
            out <- array(NA_real_, c(length(ctxs), horizon, 9))
            for (i in seq_along(ctxs)) {
                ctx <- ctxs[[i]]
                x <- ctx$y
                tab <- build_lag_table(
                    x,
                    do = ctx$do,
                    dosat = ctx$dosat,
                    temp = ctx$temp,
                    airtemp = ctx$airtemp
                )
                if (is.null(tab)) {
                    out[i, , ] <- x[length(x)]
                    next
                }
                step_fn <- fit_and_predict(tab)
                paths <- simulate_paths(
                    x,
                    step_fn,
                    horizon,
                    199,
                    tail(ctx$do, 1),
                    tail(ctx$dosat, 1),
                    tail(ctx$temp, 1),
                    tail(ctx$airtemp, 1)
                )
                out[i, , ] <- t(apply(paths, 2, quantile, probs = LEVELS, names = FALSE))
            }
            out
        }
    )
}

backend_rf <- function() {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("install.packages('ranger')")
    seed <- as.integer(Sys.getenv("SEED", "42"))
    backend_lagmodel("rf", function(tab) {
        fit <- ranger::ranger(
            y ~ .,
            data = tab,
            quantreg = TRUE,
            num.trees = 500,
            seed = seed
        )
        function(newdata) {
            predict(fit, data = newdata, type = "quantiles", quantiles = LEVELS)$predictions
        }
    })
}
