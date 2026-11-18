library(readr)
library(dplyr)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)

chla <- read_csv("chla_daily.csv", show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(Chla_ugL_mean)) %>%
  arrange(date)

n <- nrow(chla)
test_start_idx <- floor(n * 0.7) + 1
chla_test <- chla[test_start_idx:n, ] %>%
  mutate(origin = row_number())

cat("Test set:", nrow(chla_test), "rows\n")
cat("Test dates:", as.character(min(chla_test$date)), "to", as.character(max(chla_test$date)), "\n")

model_colors <- c(
  "ens3" = "#e6127c",
  "arima" = "#ff009d",
  "rf" = "#ff00a2",
  "timesfm" = "#f60000",
  "chronos" = "#ff4bf0"
)

model_names <- c(
  "ens3" = "ens3 (Ensemble)",
  "arima" = "ARIMA",
  "rf" = "Random Forest",
  "timesfm" = "TimesFM",
  "chronos" = "Chronos"
)

models <- c("ens3", "arima", "chronos", "rf", "timesfm")

for (m in models) {
  scores <- read_csv(paste0("results_all5/scores_", m, "_test_raw.csv"),
                     show_col_types = FALSE)
  min_h <- min(scores$h)

  df <- scores %>%
    filter(h == min_h) %>%
    mutate(origin = as.integer(origin)) %>%
    left_join(chla_test %>% select(origin, date), by = "origin")

  cat(m, "- joined rows:", nrow(df), "NA dates:", sum(is.na(df$date)), "\n")

  if (sum(!is.na(df$date)) == 0) next

  p <- ggplot(df %>% filter(!is.na(date)), aes(x = date)) +
    geom_line(aes(y = y), color = "gray40", linewidth = 0.8, alpha = 0.8) +
    geom_line(aes(y = med), color = model_colors[m], linewidth = 0.8) +
    labs(
      title = paste0(model_names[m], " — Actual vs Predicted (Test Set, h=", min_h, ")"),
      x = "Date", y = "Chlorophyll-a (μg/L)",
      caption = "Gray = Observed. Colored = Model median forecast."
    ) +
    theme_minimal(base_size = 13)

  ggsave(paste0("plots/test_", m, "_actual_vs_pred.png"),
         p, width = 10, height = 5, dpi = 300)
  cat("Saved test_", m, "_actual_vs_pred.png\n")
}

cat("\nDone!\n")
