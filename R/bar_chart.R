library(ggplot2)

df <- data.frame(
    Model = c("Persistence", "ARIMA", "Random Forest", "TimesFM", "Chronos", "Ensemble (ens3)"),
    CRPS = c(10.271, 8.941, 10.936, 7.942, 8.679, 7.853)
)


df$Model <- factor(df$Model, levels = df$Model[order(df$CRPS, decreasing = TRUE)])

ggplot(df, aes(x = CRPS, y = Model, fill = Model)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = round(CRPS, 3)),
        hjust = -0.2, size = 4, fontface = "bold"
    ) +
    scale_fill_manual(values = c(
        "Persistence"     = "#b0b0b0",
        "ARIMA"           = "#5B8DB8",
        "Random Forest"   = "#5B8DB8",
        "TimesFM"         = "#1B4F8A",
        "Chronos"         = "#1B4F8A",
        "Ensemble (ens3)" = "#0D2E5A"
    )) +
    labs(
        title = "Model Comparison: CRPS (Test Set, Calibrated)",
        subtitle = "Lower CRPS = Better Performance",
        x = "CRPS",
        y = NULL
    ) +
    xlim(0, 13) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray50"),
        axis.text.y = element_text(face = "bold")
    ) +
    geom_vline(
        xintercept = 10.271, linetype = "dashed",
        color = "gray40", linewidth = 0.8
    ) +
    annotate("text",
        x = 9.5, y = 0.5,
        label = "Baseline", color = "gray40", size = 3.5, hjust = 1
    )

ggsave("crps_comparison.png", width = 9, height = 7, dpi = 300)
