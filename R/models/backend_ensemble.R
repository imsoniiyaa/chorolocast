ens <- {
    w <- as.numeric(Sys.getenv("ENS_W", "0.6"))

    backend_ensemble(
        backend_timesfm(ctx_len),
        backend_rf(),
        w
    )
}

backend_ensemble_h <- function(b1, b2, w0 = 1, w1 = 0.75) {
    list(
        name = "ens_h",
        predict = function(ctxs, horizon) {
            Q1 <- b1$predict(ctxs, horizon)
            Q2 <- b2$predict(ctxs, horizon)
            w <- w0 + (w1 - w0) * pmin(1, (seq_len(horizon) - 1) / (HORIZON - 1))
            for (h in seq_len(horizon)) {
                Q1[, h, ] <- w[h] * Q1[, h, ] + (1 - w[h]) * Q2[, h, ]
            }
            Q1
        }
    )
}
backend_ens3 <- function(b_arima, b_rf, b_tfm, a0 = 1, a1 = 0.75, v0 = 0.75, v1 = 0.5) {
    list(name = "ens3", predict = function(ctxs, horizon) {
        Qa <- b_arima$predict(ctxs, horizon)
        Qr <- b_rf$predict(ctxs, horizon)
        Qt <- b_tfm$predict(ctxs, horizon)
        s <- pmin(1, (seq_len(horizon) - 1) / (HORIZON - 1))
        a <- a0 + (a1 - a0) * s
        v <- v0 + (v1 - v0) * s
        out <- Qt
        for (h in seq_len(horizon)) {
            Qe <- a[h] * Qa[, h, ] + (1 - a[h]) * Qr[, h, ]
            out[, h, ] <- v[h] * Qt[, h, ] + (1 - v[h]) * Qe
        }
        out
    })
}
