# finish a time series by filling in missing days and averaging duplicates
finish_series <- function(d, y) {
    m <- tapply(y, d, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
    days <- seq(min(as.Date(names(m))), max(as.Date(names(m))), by = "day")
    data.frame(date = days, y = as.numeric(m[as.character(days)]))
}
# load a time series from a CSV file or from VERA
load_series <- function(path, date_col = NULL, target = VARIABLE) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (is.null(date_col)) {
        date_col <- intersect(c("Date", "date", "datetime", "DateTime", "time", "timestamp"), names(df))[1]
    }
    if (is.na(date_col)) stop("no date column found in: ", paste(names(df), collapse = ", "), " (use --date-col)")
    finish_series(as.Date(substr(as.character(df[[date_col]]), 1, 10)), as.numeric(df[[target]]))
}

# converting to and from log1p space, with a floor at 0
fwd <- function(x, log) {
    x <- pmax(x, 0)
    if (log) log1p(x) else x
}
inv <- function(x, log) pmax(if (log) expm1(x) else x, 0)

# make a context vector for a given origin date, or return NULL if not enough data
make_context <- function(series, origin, ctx_len, log, max_stale = 7, min_obs = 30) {
    i <- match(origin, series$date)
    if (is.na(i)) {
        return(NULL)
    }
    s <- series$y[max(1, i - ctx_len + 1):i]
    if (sum(!is.na(s)) < min_obs) {
        return(NULL)
    }
    if (length(s) - max(which(!is.na(s))) > max_stale) {
        return(NULL)
    }
    s <- approx(seq_along(s), s, xout = seq_along(s), rule = 2)$y
    fwd(s, log)
}