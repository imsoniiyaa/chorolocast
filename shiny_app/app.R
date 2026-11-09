library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(plotly)

base_url <- "https://raw.githubusercontent.com/imsoniiyaa/chorolocast/main/"

# 색상 팔레트
colors <- list(
  primary = "#2D6A4F",
  secondary = "#52B788",
  accent = "#95D5B2",
  light = "#D8F3DC",
  bg = "#F8FFF9",
  gray = "#B0B0B0",
  ml = "#40916C",
  foundation = "#1B4332"
)

ui <- fluidPage(
  tags$head(
    tags$style(HTML(sprintf(
      "
      body { background-color: %s; font-family: 'Helvetica Neue', sans-serif; }
      .title-panel {
        background: linear-gradient(135deg, %s, %s);
        color: white; padding: 25px 30px; border-radius: 12px;
        margin-bottom: 20px;
      }
      .card {
        background: white; border-radius: 12px;
        padding: 20px; margin-bottom: 20px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      }
      .stat-box {
        background: linear-gradient(135deg, %s, %s);
        color: white; border-radius: 10px;
        padding: 15px; text-align: center;
        margin: 5px;
      }
      .stat-box h3 { font-size: 28px; margin: 0; font-weight: bold; }
      .stat-box p { margin: 0; font-size: 13px; opacity: 0.9; }
      h3 { color: %s; }
    ", colors$bg, colors$primary, colors$secondary,
      colors$secondary, colors$primary, colors$primary
    )))
  ),

  # 타이틀
  div(
    class = "title-panel",
    h2("🌿 Chlorocast Dashboard", style = "margin:0; font-size:28px;"),
    p(paste("Chlorophyll-a Forecasting | Last updated:", Sys.Date()),
      style = "margin:5px 0 0 0; opacity:0.85;"
    )
  ),

  # 상단 stat boxes
  fluidRow(
    column(3, div(class = "stat-box", h3("7.853"), p("Best CRPS (ens3)"))),
    column(3, div(class = "stat-box", h3("2,455"), p("Daily Records"))),
    column(3, div(class = "stat-box", h3("35"), p("Days Forecast Horizon"))),
    column(3, div(class = "stat-box", h3("6"), p("Models Evaluated")))
  ),

  # CRPS + Horizon
  fluidRow(
    column(
      6,
      div(
        class = "card",
        h3("Model Comparison: CRPS"),
        p("Lower CRPS = Better Performance", style = "color:gray; margin-top:-10px;"),
        plotlyOutput("crps_plot", height = "350px")
      )
    ),
    column(
      6,
      div(
        class = "card",
        h3("CRPS by Forecast Horizon"),
        p("Performance across 1-35 day horizons", style = "color:gray; margin-top:-10px;"),
        plotlyOutput("horizon_plot", height = "350px")
      )
    )
  ),

  # Forecast plot
  fluidRow(
    column(
      12,
      div(
        class = "card",
        h3("Ensemble Forecast (ens3)"),
        p("Shaded area = 90% prediction interval", style = "color:gray; margin-top:-10px;"),
        plotlyOutput("forecast_plot", height = "400px")
      )
    )
  )
)

server <- function(input, output) {
  metrics <- reactive({
    read_csv(paste0(base_url, "results_all5/metrics_overall.csv"), show_col_types = FALSE)
  })
  by_h <- reactive({
    read_csv(paste0(base_url, "results_all5/metrics_by_horizon.csv"), show_col_types = FALSE)
  })
  forecast <- reactive({
    read_csv(paste0(base_url, "forecasts/chlorocast_ens3.csv"), show_col_types = FALSE)
  })

  output$crps_plot <- renderPlotly({
    df <- metrics() %>%
      filter(split == "test", variant == "calibrated") %>%
      mutate(
        model = factor(model, levels = model[order(CRPS, decreasing = TRUE)]),
        color = case_when(
          model == "ens3" ~ colors$foundation,
          model %in% c("timesfm", "chronos") ~ colors$primary,
          model == "persistence" ~ colors$gray,
          TRUE ~ colors$secondary
        )
      )

    p <- ggplot(df, aes(
      x = CRPS, y = model, fill = model,
      text = paste0("Model: ", model, "<br>CRPS: ", round(CRPS, 3))
    )) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_vline(
        xintercept = df$CRPS[df$model == "persistence"],
        linetype = "dashed", color = "gray50"
      ) +
      scale_fill_manual(values = setNames(df$color, df$model)) +
      xlim(0, 13) +
      labs(x = "CRPS", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    ggplotly(p, tooltip = "text")
  })

  output$horizon_plot <- renderPlotly({
    df <- by_h() %>%
      filter(split == "test", variant == "calibrated")

    p <- ggplot(df, aes(
      x = h, y = CRPS, color = model, group = model,
      text = paste0("Model: ", model, "<br>Horizon: ", h, " days<br>CRPS: ", round(CRPS, 3))
    )) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c(
        "persistence" = colors$gray,
        "arima" = colors$accent,
        "rf" = colors$secondary,
        "timesfm" = colors$primary,
        "chronos" = "#74C69D",
        "ens3" = colors$foundation
      )) +
      labs(x = "Forecast Horizon (days)", y = "CRPS", color = "Model") +
      theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    ggplotly(p, tooltip = "text")
  })

  output$forecast_plot <- renderPlotly({
    df <- forecast() %>%
      filter(parameter %in% c(10, 100, 190)) %>%
      mutate(datetime = as.Date(datetime)) %>%
      pivot_wider(names_from = parameter, values_from = prediction) %>%
      rename(lo = `10`, mid = `100`, hi = `190`)

    plot_ly(df, x = ~datetime) %>%
      add_ribbons(
        ymin = ~lo, ymax = ~hi,
        fillcolor = "rgba(82, 183, 136, 0.2)",
        line = list(color = "transparent"),
        name = "90% Interval"
      ) %>%
      add_lines(
        y = ~mid,
        line = list(color = colors$primary, width = 2),
        name = "Forecast (median)",
        hovertemplate = "Date: %{x}<br>Chla: %{y:.2f} μg/L<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = "Date"),
        yaxis = list(title = "Chlorophyll-a (μg/L)"),
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        legend = list(orientation = "h")
      )
  })
}

shinyApp(ui = ui, server = server)