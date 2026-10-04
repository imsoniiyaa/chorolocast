if (!exists("LEVELS")) {
    LEVELS <- seq(0.1, 0.9, by = 0.1)
}

MID <- which.min(abs(LEVELS - 0.5))

backend_log <- function(backend, back = TRUE) {
    list(
        name = paste0(backend$name, "_log"),
        predict = function(ctxs, horizon) {
            ctxs_l <- lapply(
                ctxs,
                function(x) {
                    log1p(pmax(x, 0))
                }
            )

            Q <- backend$predict(
                ctxs_l,
                horizon
            )

            if (back) {
                expm1(Q)
            } else {
                Q
            }
        }
    )
}

backend_ensemble <- function(b1, b2, w) {
    list(
        name = "ens",
        predict = function(ctxs, horizon) {
            w * b1$predict(ctxs, horizon) +
                (1 - w) * b2$predict(ctxs, horizon)
        }
    )
}

conformal_quantiles <- function(
  val_Q,
  val_y,
  test_Q,
  pool = 3
) {
    H <- dim(test_Q)[2]

    out <- test_Q

    for (h in seq_len(H)) {
        hs <- max(1, h - pool):min(H, h + pool)

        res <- as.vector(
            val_y[, hs] -
                val_Q[, hs, MID]
        )

        res <- res[is.finite(res)]

        adj <- quantile(
            res,
            probs = LEVELS,
            names = FALSE,
            type = 1
        )

        out[, h, ] <-
            test_Q[, h, MID] +
            matrix(
                adj,
                nrow = dim(test_Q)[1],
                ncol = length(LEVELS),
                byrow = TRUE
            )
    }

    out
}

pinball <- function(Q, y) {
    mean(
        sapply(
            seq_along(LEVELS),
            function(k) {
                e <- y - Q[, , k]

                mean(
                    pmax(
                        LEVELS[k] * e,
                        (LEVELS[k] - 1) * e
                    ),
                    na.rm = TRUE
                )
            }
        )
    )
}

blend <- function(Qa, Qb, w) {
    w * Qa + (1 - w) * Qb
}

pick_weight <- function(
  val_Qa,
  val_Qb,
  val_y,
  ws = seq(0, 1, 0.01)
) {
    scores <- sapply(
        ws,
        function(w) {
            pinball(
                blend(val_Qa, val_Qb, w),
                val_y
            )
        }
    )

    ws[which.min(scores)]
}