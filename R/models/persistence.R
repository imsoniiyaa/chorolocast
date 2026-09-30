# persistence model: predict that the future will be the same as the last observed value, 
# plus a quantile of the recent changes
backend_persistence <- function() {
  list(name = "persistence", predict = function(ctxs, horizon) {
    out <- array(0, c(length(ctxs), horizon, 9))
    for (i in seq_along(ctxs)) {
      x <- ctxs[[i]]; n <- length(x)
      for (h in seq_len(horizon)) {
        d <- if (n - h >= 5) x[(h + 1):n] - x[1:(n - h)] else NULL
        out[i, h, ] <- x[n] + if (is.null(d)) 0 else quantile(d, LEVELS, names = FALSE)
      }
    }
    out
  })
}