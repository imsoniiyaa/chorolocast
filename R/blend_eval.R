source("R/utils/ensemble.R")
source("R/utils/python_backend.R")
source("R/utils/data.R")
source("R/utils/calibration.R")
source("R/utils/metrics.R")
source("R/utils/forecast.R")
source("R/models/backend_ensemble.R")
source("R/models/persistence.R")
source("R/models/arima.R")
source("R/models/rf.R")

#Arima + RF ensemble

HORIZON <- 35L
LEVELS <- seq(0.1, 0.9, by = 0.1)
Z90 <- qnorm(0.9)
VARIABLE <- "Chla_ugL_mean"

series <- load_series("chla_daily.csv", NULL, VARIABLE)
by <- "3 days"
origins <- list(
    val = seq(as.Date("2025-01-01"), as.Date("2025-12-31"), by = by),
    test = seq(as.Date("2026-01-01"), as.Date("2026-09-20"), by = by)
)

b_a <- backend_arima()
b_r <- backend_rf()
R <- lapply(origins, function(o) {
    list(
        a = collect(series, b_a, o, 512, TRUE),
        r = collect(series, b_r, o, 512, TRUE)
    )
})

blend_h <- function(Qa, Qr, w) {
    out <- Qa
    for (h in seq_len(dim(Qa)[2])) out[, h, ] <- w[h] * Qa[, h, ] + (1 - w[h]) * Qr[, h, ]
    out
}
wvec <- function(w0, w1) w0 + (w1 - w0) * (seq_len(HORIZON) - 1) / (HORIZON - 1)
crps_of <- function(sp, w0, w1) {
    Q <- blend_h(R[[sp]]$a$Q, R[[sp]]$r$Q, wvec(w0, w1))
    mean(score_forecasts(Q, R[[sp]]$a$Y, TRUE, NULL, 199)$crps)
}

g <- expand.grid(w0 = seq(0, 1, 0.25), w1 = seq(0, 1, 0.25))
g$val <- mapply(function(a, b) crps_of("val", a, b), g$w0, g$w1)
g <- g[order(g$val), ]
cat("\n=== val CRPS 상위 8개 (w = arima 비중, w0: h=1, w1: h=35) ===\n")
print(head(g, 8), row.names = FALSE)

best <- g[1, ]
cat("\n=== test CRPS ===\n")
res <- data.frame(
    name = c("arima only", "rf only", "best blend"),
    w0 = c(1, 0, best$w0), w1 = c(1, 0, best$w1)
)
res$val <- mapply(function(a, b) crps_of("val", a, b), res$w0, res$w1)
res$test <- mapply(function(a, b) crps_of("test", a, b), res$w0, res$w1)
print(res, row.names = FALSE)