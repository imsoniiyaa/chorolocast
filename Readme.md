# Chlorocast

**Chlorophyll-a Forecasting for Falling Creek Reservoir**  
Virginia Tech CMDA Capstone Project | Center for Ecosystem Forecasting (CEF)

[![Dashboard](https://img.shields.io/badge/Dashboard-Live-brightgreen)](https://chorolocast.shinyapps.io/shiny_app/)

---

## Overview

Chlorocast is an ensemble forecasting system for chlorophyll-a (Chla) concentrations in Falling Creek Reservoir, Virginia. Harmful algal blooms driven by high Chla levels pose risks to drinking water quality. This project produces daily probabilistic forecasts up to 35 days ahead, submitted to the **VERA Forecasting Challenge**.

---

## Dashboard

🔗 [https://chorolocast.shinyapps.io/shiny_app/](https://chorolocast.shinyapps.io/shiny_app/)

- Live CRPS model comparison
- Forecast horizon performance (1–35 days)
- Daily ensemble forecast with 90% prediction interval

---

## Models

| Model         | Description                                    | Test CRPS |
| ------------- | ---------------------------------------------- | --------- |
| Persistence   | Baseline: last observed value                  | 10.271    |
| ARIMA         | Statistical time series                        | 8.941     |
| Random Forest | ML with lagged features                        | 10.936    |
| TimesFM       | Google foundation model (zero-shot)            | 7.942     |
| Chronos       | Amazon foundation model (zero-shot)            | 8.679     |
| **ens3**      | **ARIMA + RF + TimesFM ensemble (calibrated)** | **6.383** |

**ens3** is the best-performing model, combining classical and foundation model forecasts with post-hoc calibration (coverage95 = 94.4%).

---

## Methods

### Data

- Source: VERA Forecast — daily in-situ sensor data
- Target: `Chla_ugL_mean` (chlorophyll-a, μg/L)
- Period: 2013–present, Falling Creek Reservoir
- Features: water temp, DO, fDOM, turbidity, air temp, PAR, wind, rain, inflow

### Evaluation

- Metric: **CRPS** (Continuous Ranked Probability Score) — lower is better
- Train/Val/Test split with stride-3 rolling window
- Calibration via isotonic regression on validation set

### Automation

- GitHub Actions runs daily at midnight UTC
- Fetches latest VERA data → runs ens3 forecast → pushes to GitHub → Shiny app auto-updates

---

## Repository Structure

chlorocast/
├── R/
│ ├── fm.R # Main R forecasting script (evaluate + forecast)
│ ├── update_data.R # Fetches latest VERA data
│ └── utils/
│ └── data.R # Data loading utilities
├── py/
│ ├── ens4.py # Python ensemble (ChloroFM + TimesFM)
│ ├── chlorofm.py # ChloroFM model
│ ├── fm_backend.py # TimesFM/Chronos backend
│ └── ...
├── config/
│ └── chlorofm.json # Model configuration
├── shiny_app/
│ └── app.R # Shiny dashboard
├── forecasts/
│ └── chlorocast_ens3.csv # Latest forecast output
├── results_all5/
│ ├── metrics_overall.csv # Model evaluation results
│ └── metrics_by_horizon.csv
└── .github/workflows/
└── forecast.yml # Daily automation

---

## Quickstart

### R (ens3)

```bash
Rscript R/fm.R forecast \
  --data chla_daily.csv \
  --model ens3 \
  --out forecasts
```

### Python (ens4)

```bash
# Requires Python 3.12, timesfm==2.0.2
python py/ens4.py forecast \
  --data chla_daily.csv \
  --out forecasts
```

---

## Results

| Split             | Model    | CRPS      | Coverage95 |
| ----------------- | -------- | --------- | ---------- |
| Val               | ens3     | 5.807     | 90.3%      |
| Test (raw)        | ens3     | 7.782     | 89.7%      |
| Test (calibrated) | **ens3** | **6.383** | **94.4%**  |

---

## References

- [VERA Forecasting Challenge](https://www.ltreb-reservoirs.org/vera4cast/)
- [TimesFM (Google)](https://github.com/google-research/timesfm)
- [Chronos (Amazon)](https://github.com/amazon-science/chronos-forecasting)
- [Center for Ecosystem Forecasting, Virginia Tech](https://www.ltreb-reservoirs.org/)
