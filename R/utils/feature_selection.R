select_covars <- function(df) {
    df[, c("lake_DO_mgL_mean", "lake_Temp_C_mean", "lake_DOsat_percent_mean", "met_AirTemp_C_mean"), drop = FALSE]
}
