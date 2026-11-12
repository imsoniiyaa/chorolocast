"""ChloroFM: TimesFM plus causal environmental/seasonal analog trajectories."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fm import HORIZON, inv, make_context, quantiles_to_samples, to_vera  # noqa: E402
from fm_backend import TimesFM25  # noqa: E402

TARGET = "Chla_ugL_mean"
ENV_COLUMNS = [
    "lake_Temp_C_mean", "lake_SpCond_uScm_mean", "lake_fDOM_QSU_mean",
    "lake_DOsat_percent_mean", "inflow_Flow_cms_mean", "inflow_Temp_C_mean",
    "met_PAR_umolm2s_mean", "met_AirTemp_C_mean", "met_RH_percent_mean",
    "met_WindSpeed_ms_mean", "met_Rain_mm_sum",
]


def load_daily(path: str | Path, date_col: str = "date") -> pd.DataFrame:
    d = pd.read_csv(path, parse_dates=[date_col]).set_index(date_col).sort_index()
    d.index = d.index.tz_convert(None)
    d.index = d.index.normalize()  # 날짜만 (00:00:00으로 정규화)
    d = d.groupby(level=0).mean(numeric_only=True).asfreq("D")
    return d


def make_features(d: pd.DataFrame) -> pd.DataFrame:
    """Build origin features causally: missing values use past observations only."""
    ch = d[TARGET].ffill()
    f = pd.DataFrame(index=d.index)
    doy = d.index.dayofyear.to_numpy()
    f["sin"] = np.sin(2 * np.pi * doy / 365.25)
    f["cos"] = np.cos(2 * np.pi * doy / 365.25)
    for window in (1, 7, 14, 30, 60):
        f[f"ch_mean{window}"] = ch.rolling(window, min_periods=1).mean()
    f["ch_last"] = ch
    f["ch_slope7"] = ch - ch.shift(7)
    for col in ENV_COLUMNS:
        f[col] = d[col].ffill().rolling(7, min_periods=1).mean()
    return f


def feature_columns(group: str, all_columns: list[str]) -> list[str]:
    if group == "season":
        return ["sin", "cos"]
    if group == "season_chla":
        return ["sin", "cos", "ch_last", "ch_mean7", "ch_mean30", "ch_slope7"]
    if group == "season_env":
        return all_columns
    raise ValueError(f"unknown analog group: {group}")


def analog_samples(d: pd.DataFrame, f: pd.DataFrame, origin: pd.Timestamp,
                   group: str, neighbors: int, recency: float,
                   members: int = 199, horizon: int = HORIZON) -> np.ndarray:
    """Return members x horizon trajectories whose outcomes were known at origin."""
    cutoff = origin - pd.Timedelta(days=horizon)
    candidates = f.index[(f.index >= f.index.min() + pd.Timedelta(days=60)) & (f.index <= cutoff)]
    if len(candidates) < 2:
        raise ValueError(f"not enough historical analogs before {origin.date()}")
    cols = feature_columns(group, list(f.columns))
    train = f.loc[candidates, cols]
    scale = train.std().replace(0, 1).fillna(1)
    distance = np.sqrt((((train - f.loc[origin, cols]) / scale) ** 2).mean(axis=1))
    if recency:
        age_years = (origin - candidates).days / 365.25
        distance *= 1 + recency * np.maximum(age_years - 1, 0)
    nearest = distance.nsmallest(min(neighbors, len(distance))).index
    target = d[TARGET].ffill().bfill()
    trajectories = np.stack([
        [target.loc[t + pd.Timedelta(days=h)] for h in range(1, horizon + 1)]
        for t in nearest
    ])
    pick = np.floor(np.linspace(0, len(trajectories), members, endpoint=False)).astype(int)
    return trajectories[pick]


def crps_ensemble(samples: np.ndarray, outcomes: np.ndarray) -> np.ndarray:
    """CRPS for origins x members x horizons and origins x horizons."""
    ordered = np.sort(samples, axis=1)
    n = ordered.shape[1]
    first = np.nanmean(np.abs(ordered - outcomes[:, None, :]), axis=1)
    weights = 2 * np.arange(1, n + 1) - n - 1
    second = np.sum(ordered * weights[None, :, None], axis=1) / n**2
    result = first - second
    result[np.isnan(outcomes)] = np.nan
    return result


def foundation_samples(backend: TimesFM25, series: pd.Series, origins: pd.DatetimeIndex,
                       members: int, context_length: int) -> np.ndarray:
    contexts = [make_context(series, o, context_length, True) for o in origins]
    quantiles = backend.predict(contexts, HORIZON)
    return np.stack([np.clip(inv(quantiles_to_samples(q, members), True), 0, None) for q in quantiles])


def hybrid_samples(d: pd.DataFrame, origins: pd.DatetimeIndex, foundation: np.ndarray,
                   config: dict, members: int) -> np.ndarray:
    f = make_features(d)
    result = np.empty_like(foundation)
    cache = {}
    for band in config["bands"]:
        key = (band["analog"], band["neighbors"], band["recency"])
        if key not in cache:
            cache[key] = np.stack([
                analog_samples(d, f, o, *key, members=members) for o in origins
            ])
        hs = slice(band["start"] - 1, band["end"])
        analog = np.sort(cache[key], axis=1)
        weight = band["timesfm_weight"]
        result[:, :, hs] = weight * np.sort(foundation, axis=1)[:, :, hs] + (1 - weight) * analog[:, :, hs]
    return np.clip(result, 0, None)


def outcomes(d: pd.DataFrame, origins: pd.DatetimeIndex) -> np.ndarray:
    result = np.full((len(origins), HORIZON), np.nan)
    for i, origin in enumerate(origins):
        result[i] = d[TARGET].reindex(pd.date_range(origin + pd.Timedelta(days=1), periods=HORIZON))
    return result


def load_config(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def cmd_evaluate(a):
    d = load_daily(a.data)
    config = load_config(a.config)
    origins = pd.date_range(a.start, a.end, freq=f"{a.stride}D")
    y = outcomes(d, origins)
    backend = TimesFM25(config["context_length"])
    base = foundation_samples(backend, d[TARGET], origins, a.members, config["context_length"])
    github_base = foundation_samples(backend, d[TARGET], origins, a.members, 512)
    hybrid = hybrid_samples(d, origins, base, config, a.members)
    rows = []
    for name, samples in [("timesfm_log512_github", github_base),
                          ("timesfm_log1024", base), ("chlorofm", hybrid)]:
        scores = crps_ensemble(samples, y)
        for h in range(HORIZON):
            rows.append({"model": name, "h": h + 1, "n": int(np.isfinite(scores[:, h]).sum()),
                         "CRPS": float(np.nanmean(scores[:, h]))})
    by_h = pd.DataFrame(rows)
    overall = by_h.groupby("model", as_index=False).apply(
        lambda x: pd.Series({"CRPS": np.average(x.CRPS, weights=x.n), "n": int(x.n.sum())}),
        include_groups=False,
    ).reset_index(drop=True)
    github_crps = float(overall.loc[overall.model == "timesfm_log512_github", "CRPS"].iloc[0])
    stronger_crps = float(overall.loc[overall.model == "timesfm_log1024", "CRPS"].iloc[0])
    overall["CRPSS_vs_GitHub_TimesFM512"] = 1 - overall.CRPS / github_crps
    overall["CRPSS_vs_TimesFM1024"] = 1 - overall.CRPS / stronger_crps
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    by_h.to_csv(out / "metrics_by_horizon.csv", index=False)
    overall.to_csv(out / "metrics_overall.csv", index=False)
    print(overall.round(4).to_string(index=False))
    print(f"saved evaluation to {out}")


def cmd_forecast(a):
    d = load_daily(a.data)
    config = load_config(a.config)
    ref = pd.Timestamp(a.reference_date) if a.reference_date else d[TARGET].last_valid_index()
    if ref not in d.index:
        raise ValueError(f"reference date {ref.date()} is outside the daily data")
    backend = TimesFM25(config["context_length"])
    origins = pd.DatetimeIndex([ref])
    base = foundation_samples(backend, d[TARGET], origins, a.members, config["context_length"])
    samples = hybrid_samples(d, origins, base, config, a.members)[0]
    forecast = to_vera(samples, ref, a.model_id or config["model_id"], a.site, a.depth, TARGET)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"daily-{ref.date()}-{forecast.model_id.iat[0]}.csv.gz"
    forecast.to_csv(path, index=False)
    print(f"wrote {path} ({len(forecast)} rows)")


def main():
    default_config = str(HERE.parent / "config" / "chlorofm.json")
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    e = sub.add_parser("evaluate")
    e.add_argument("--data", required=True)
    e.add_argument("--config", default=default_config)
    e.add_argument("--start", default="2026-01-01")
    e.add_argument("--end", default="2026-08-16")
    e.add_argument("--stride", type=int, default=7)
    e.add_argument("--members", type=int, default=199)
    e.add_argument("--out", default="results_chlorofm")
    e.set_defaults(fn=cmd_evaluate)
    f = sub.add_parser("forecast")
    f.add_argument("--data", required=True)
    f.add_argument("--config", default=default_config)
    f.add_argument("--reference-date")
    f.add_argument("--members", type=int, default=199)
    f.add_argument("--model-id")
    f.add_argument("--site", default="fcre")
    f.add_argument("--depth", type=float, default=1.6)
    f.add_argument("--out", default="forecasts")
    f.set_defaults(fn=cmd_forecast)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
