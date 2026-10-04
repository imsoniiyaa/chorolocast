# finish a time series by filling in missing days and averaging duplicates
finish_series <- function(df, date_col, target) {
    df$date <- as.Date(
        substr(as.character(df[[date_col]]), 1, 10)
    )

    keep <- c(
        "date",
        target,
        "lake_DO_mgL_mean",
        "lake_DOsat_percent_mean",
        "lake_Temp_C_mean",
        "met_AirTemp_C_mean"
    )

    df <- df[, keep]

    names(df)[names(df) == target] <- "y"

    df
}
# load a time series from a CSV file or from VERA
load_series <- function(path,
                        date_col = NULL,
                        target = VARIABLE) {
    df <- read.csv(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (is.null(date_col)) {
        date_col <- intersect(
            c(
                "Date",
                "date",
                "datetime",
                "DateTime",
                "time",
                "timestamp"
            ),
            names(df)
        )[1]
    }

    finish_series(
        df,
        date_col,
        target
    )
}

# converting to and from log1p space, with a floor at 0
fwd <- function(x, log) {
    x <- pmax(x, 0)
    if (log) log1p(x) else x
}
inv <- function(x, log) pmax(if (log) expm1(x) else x, 0)

# make a context vector for a given origin date, or return NULL if not enough data

make_context <- function(series, origin, ctx_len, log,
                         max_stale = 7, min_obs = 30) {
    i <- match(origin, series$date)

    if (is.na(i)) {
        return(NULL)
    }

    win <- max(1, i - ctx_len + 1):i

    s <- series$y[win]

    if (sum(!is.na(s)) < min_obs) {
        return(NULL)
    }

    if (length(s) - max(which(!is.na(s))) > max_stale) {
        return(NULL)
    }

    s <- approx(
        seq_along(s),
        s,
        xout = seq_along(s),
        rule = 2
    )$y

    list(
        y = fwd(s, log),
        do = series$lake_DO_mgL_mean[win],
        dosat = series$lake_DOsat_percent_mean[win],
        temp = series$lake_Temp_C_mean[win],
        airtemp = series$met_AirTemp_C_mean[win]
    )
}
