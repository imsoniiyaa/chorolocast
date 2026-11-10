library(readr)
library(dplyr)
library(tidyr)

url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"
url_inflow <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-inflow-targets.csv.gz"

targets <- read_csv(url, show_col_types = FALSE)
inflow_targets <- read_csv(url_inflow, show_col_types = FALSE)
targets <- dplyr::bind_rows(targets, inflow_targets)

chla <- targets %>%
    filter(site_id == "fcre") %>%
    pivot_wider(
        id_cols = c(datetime, site_id, depth_m),
        names_from = variable,
        values_from = observation,
        values_fn = mean
    ) %>%
    filter(depth_m == 1.6 | is.na(depth_m)) %>%
    arrange(datetime) %>%
    rename(date = datetime)

# 없는 컬럼은 NA로 채우기
required_cols <- c(
    "date", "Chla_ugL_mean", "lake_DO_mgL_mean",
    "lake_DOsat_percent_mean", "lake_Temp_C_mean",
    "met_AirTemp_C_mean", "lake_SpCond_uScm_mean",
    "lake_fDOM_QSU_mean", "lake_Turbidity_FNU_mean",
    "inflow_Flow_cms_mean", "inflow_Temp_C_mean",
    "met_PAR_umolm2s_mean", "met_BP_kPa_mean",
    "met_RH_percent_mean", "met_WindSpeed_ms_mean",
    "met_WindDir_degrees_mean", "met_Rain_mm_sum"
)

for (col in required_cols) {
    if (!col %in% names(chla)) chla[[col]] <- NA
}

chla <- chla %>% select(all_of(required_cols))

write_csv(chla, "chla_daily.csv")
cat("Updated chla_daily.csv:", nrow(chla), "rows\n")
cat("Columns:", paste(names(chla), collapse = ", "), "\n")