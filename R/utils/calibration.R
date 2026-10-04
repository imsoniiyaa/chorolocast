# convert a matrix of quantiles to samples,
# using linear interpolation for the inner quantiles and normal approximation for the tails
quantiles_to_samples <- function(q, n) {
    q <- matrix(q, ncol = 9)
    q <- t(apply(q, 1, sort))
    H <- nrow(q)
    u <- (seq_len(n) - 0.5) / n
    inner <- u >= 0.1 & u <= 0.9
    lo <- u < 0.1
    hi <- u > 0.9
    s_lo <- pmax((q[, 5] - q[, 1]) / Z90, 1e-6)
    s_hi <- pmax((q[, 9] - q[, 5]) / Z90, 1e-6)
    out <- matrix(0, n, H)
    for (h in seq_len(H)) {
        out[inner, h] <- approx(LEVELS, q[h, ], xout = u[inner], rule = 2)$y
        out[lo, h] <- q[h, 1] + s_lo[h] * (qnorm(u[lo]) - qnorm(0.1))
        out[hi, h] <- q[h, 9] + s_hi[h] * (qnorm(u[hi]) - qnorm(0.9))
    }
    out
}

# compute the CRPS for an ensemble of samples and a vector of observations
apply_scale <- function(samples, scale) {
    med <- apply(samples, 2, median)
    scale <- scale[seq_len(ncol(samples))]
    sweep(sweep(sweep(samples, 2, med, "-"), 2, scale, "*"), 2, med, "+")
}

fit_scale <- function(Q, Y, n = 199, cover = 0.95, smooth = 5,
                      mode = Sys.getenv("CAL_MODE", "global")) {
    O <- dim(Q)[1]
    H <- dim(Q)[2]
    scores <- matrix(NA_real_, O, H)
    for (o in seq_len(O)) {
        s <- quantiles_to_samples(matrix(Q[o, , ], nrow = H), n)
        med <- apply(s, 2, median)
        ci <- apply(s, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
        scores[o, ] <- abs(Y[o, ] - med) / pmax((ci[2, ] - ci[1, ]) / 2, 1e-6)
    }

    clamp <- function(v) pmin(pmax(v, 0.25), 6)
    all_v <- scores[!is.na(scores)]
    if (length(all_v) < 10) {
        return(rep(1, H))
    }

    if (mode == "global") {
        return(rep(clamp(unname(quantile(all_v, cover))), H))
    }

    if (mode == "linear") {
        cs <- vapply(seq_len(H), function(h) {
            v <- scores[, max(1, h - 3):min(H, h + 3)]
            v <- v[!is.na(v)]
            if (length(v) >= 10) unname(quantile(v, cover)) else NA_real_
        }, numeric(1))
        ok <- !is.na(cs)
        if (sum(ok) < 3) {
            return(rep(clamp(unname(quantile(all_v, cover))), H))
        }
        fit <- lm(cs[ok] ~ which(ok))
        return(clamp(unname(coef(fit)[1] + coef(fit)[2] * seq_len(H))))
    }

    # per_h
    cs <- vapply(seq_len(H), function(h) {
        v <- scores[, h]
        v <- v[!is.na(v)]
        if (length(v) >= 10) unname(quantile(v, cover)) else NA_real_
    }, numeric(1))
    ok <- !is.na(cs)
    if (!any(ok)) {
        return(rep(1, H))
    }
    if (sum(ok) == 1) cs <- rep(cs[ok], H) else cs <- approx(which(ok), cs[ok], xout = seq_len(H), rule = 2)$y
    half <- smooth %/% 2
    cs <- vapply(seq_len(H), function(h) median(cs[max(1, h - half):min(H, h + half)]), numeric(1))
    clamp(cs)
}
