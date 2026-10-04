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
source("R/models/timesfm.R")

# TimesFM + RF + Arima ensemble

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
b_t <- backend_timesfm(512)
R <- lapply(origins, function(o) {
    list(
        a = collect(series, b_a, o, 512, TRUE),
        r = collect(series, b_r, o, 512, TRUE),
        t = collect(series, b_t, o, 512, TRUE)
    )
})

lin <- function(a, b) a + (b - a) * (seq_len(HORIZON) - 1) / (HORIZON - 1)
mix <- function(Q1, Q2, w) {
    out <- Q1
    for (h in seq_len(dim(Q1)[2])) out[, h, ] <- w[h] * Q1[, h, ] + (1 - w[h]) * Q2[, h, ]
    out
}
final_Q <- function(sp, v0, v1) {
    Qe <- mix(R[[sp]]$a$Q, R[[sp]]$r$Q, lin(1, 0.75)) 
    mix(R[[sp]]$t$Q, Qe, lin(v0, v1)) 
}
crps_of <- function(sp, v0, v1) {
    mean(score_forecasts(final_Q(sp, v0, v1), R[[sp]]$t$Y, TRUE, NULL, 199)$crps)
}

g <- expand.grid(
    v0 = seq(0.70, 1.00, 0.05),
    v1 = seq(0.40, 0.80, 0.05)
)
g$val <- mapply(function(a, b) crps_of("val", a, b), g$v0, g$v1)
g <- g[order(g$val), ]
cat("\n=== val CRPS top 8 (v = timesfm, v0: h=1, v1: h=35) ===\n")
print(head(g, 8), row.names = FALSE)

best <- g[1, ]
res <- data.frame(
    name = c("ens_h only", "timesfm only", "best blend"),
    v0 = c(0, 1, best$v0), v1 = c(0, 1, best$v1)
)
res$val <- mapply(function(a, b) crps_of("val", a, b), res$v0, res$v1)
res$test <- mapply(function(a, b) crps_of("test", a, b), res$v0, res$v1)
cat("\n=== test CRPS ===\n")
print(res, row.names = FALSE)