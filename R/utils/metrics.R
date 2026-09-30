# compute the CRPS for an ensemble of samples and a vector of observations
crps_ens <- function(S, y) {
    n <- nrow(S)
    Ss <- matrix(apply(S, 2, sort), nrow = n)
    i <- seq_len(n)
    colMeans(abs(sweep(Ss, 2, y, "-"))) - colSums((2 * i - n - 1) * Ss) / n^2
}

# compute the forecast scores for a set of quantiles and observations, optionally applying a scale
score_forecasts <- function(Q, Y, log, scale = NULL, n = 199) {
    O <- dim(Q)[1]
    H <- dim(Q)[2]
    rows <- vector("list", O)
    for (o in seq_len(O)) {
        s <- quantiles_to_samples(matrix(Q[o, , ], nrow = H), n)
        if (!is.null(scale)) s <- apply_scale(s, scale)
        S <- inv(s, log)
        y <- Y[o, ]
        ci <- apply(S, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
        med <- apply(S, 2, median)
        hh <- which(!is.na(y))
        cr <- crps_ens(S, ifelse(is.na(y), 0, y))
        rows[[o]] <- data.frame(
            origin = o, h = hh, y = y[hh], med = med[hh], crps = cr[hh],
            cover95 = as.numeric(ci[1, hh] <= y[hh] & y[hh] <= ci[2, hh]),
            width95 = ci[2, hh] - ci[1, hh]
        )
    }
    d <- do.call(rbind, rows)
    d$abs_err <- abs(d$y - d$med)
    d$sq_err <- (d$y - d$med)^2
    d
}
# additional summary functions for the forecast scores, by horizon or overall
summarise_h <- function(d) {
    do.call(rbind, lapply(split(d, d$h), function(x) {
        data.frame(
            h = x$h[1], n = nrow(x), MAE = mean(x$abs_err), RMSE = sqrt(mean(x$sq_err)),
            CRPS = mean(x$crps), coverage95 = mean(x$cover95), width95 = mean(x$width95)
        )
    }))
}
# summarize the forecast scores overall, including separate metrics for bloom events
summarise_all <- function(d) {
    data.frame(
        n = nrow(d), MAE = mean(d$abs_err), RMSE = sqrt(mean(d$sq_err)), CRPS = mean(d$crps),
        coverage95 = mean(d$cover95), width95 = mean(d$width95)
    )
}
