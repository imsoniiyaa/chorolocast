#Work same as R just to verify with python
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
#Setting standard values for the model
HORIZON = 35
LEVELS = np.round(np.arange(1, 10) / 10, 1)
MED = 4
VARIABLE = "Chla_ugL_mean"
Z90 = norm.ppf(0.9)

#get the data(cleaned earlier) and load it into a pandas dataframe. 
#If the date column is not specified, it will try to find a date column in the dataframe. 
#If it cannot find one, it will raise an error. It will then group the data by date and return a series of the target variable.
def load_series(path, date_col=None, target=VARIABLE):
    df = pd.read_csv(path)
    if date_col is None:
        date_col = next(
            (c for c in df.columns if c.lower() in ("date", "datetime", "time", "timestamp")), None
        )
        if date_col is None:
            raise SystemExit(f"Could not find a date column in {list(df.columns)}; use --date-col")
    d = pd.to_datetime(df[date_col])
    if getattr(d.dt, "tz", None) is not None:
        d = d.dt.tz_convert("UTC").dt.tz_localize(None)
    df[date_col] = d.dt.normalize()
    return df.groupby(date_col)[target].mean().asfreq("D")

#change numeric values to log scale if log is True, else return the original values.
def fwd(x, log):
    x = np.clip(x, 0, None)
    return np.log1p(x) if log else x

#change log scale values to numeric values if log is True, else return the original values.
def inv(x, log):
    return np.clip(np.expm1(x) if log else x, 0, None)

#check if there are any value to use for prediction
#If there are no values, it will return None. 
#If there are values, it will check if the last value is within the max_stale days of the origin. 
def make_context(series, origin, ctx_len, log, max_stale=7, min_obs=30):
    s = series.loc[:origin].iloc[-ctx_len:]
    if s.notna().sum() < min_obs:
        return None
    last = s.last_valid_index()
    if (origin - last).days > max_stale:
        return None
    s = s.interpolate(limit_direction="both")
    return fwd(s.to_numpy(dtype=np.float64), log).astype(np.float32)

#Simple prediction model
#This model will take the last value of the context and repeat it for the horizon.
class Persistence:

    name = "persistence"

    def predict(self, ctxs, horizon):
        out = np.empty((len(ctxs), horizon, 9), dtype=np.float32)
        for i, c in enumerate(ctxs):
            for h in range(1, horizon + 1):
                d = c[h:] - c[:-h]
                out[i, h - 1] = c[-1] + (np.quantile(d, LEVELS) if len(d) >= 5 else 0.0)
        return out

#Choose which model to use for prediction in persistence, choronos or timesfm.
def get_backend(name, ctx_len):
    if name == "persistence":
        return Persistence()
    from fm_backend import load_backend 

    return load_backend(name, ctx_len)

#fill the missing values in the context with the last available value.
def quantiles_to_samples(q, n):
    q = np.sort(q, axis=1)
    H = q.shape[0]
    u = (np.arange(n) + 0.5) / n
    inner, lo, hi = (u >= 0.1) & (u <= 0.9), u < 0.1, u > 0.9
    s_lo = np.maximum((q[:, MED] - q[:, 0]) / Z90, 1e-6)
    s_hi = np.maximum((q[:, 8] - q[:, MED]) / Z90, 1e-6)
    out = np.empty((n, H))
    for h in range(H):
        out[inner, h] = np.interp(u[inner], LEVELS, q[h])
        out[lo, h] = q[h, 0] + s_lo[h] * (norm.ppf(u[lo]) - norm.ppf(0.1))
        out[hi, h] = q[h, 8] + s_hi[h] * (norm.ppf(u[hi]) - norm.ppf(0.9))
    return out

#make sure the scale of the values is correct. 
#If the scale is not provided, it will return the original values.
def apply_scale(samples, c):
    med = np.median(samples, axis=0)
    return med + (samples - med) * np.asarray(c)[None, :samples.shape[1]]

#make the scale of the values smooth. 
#It will take the median of the values in a window of size smooth.
def fit_scale(Q, Y, n=199, cover=0.95, smooth=5):
    """Per-horizon factor c_h so that median +/- c_h * raw half-width covers ~95% (split conformal)."""
    O, H, _ = Q.shape
    scores = np.full((O, H), np.nan)
    for o in range(O):
        s = quantiles_to_samples(Q[o], n)
        med = np.median(s, axis=0)
        lo, hi = np.quantile(s, [0.025, 0.975], axis=0)
        scores[o] = np.abs(Y[o] - med) / np.maximum((hi - lo) / 2, 1e-6)
    c = np.array(
        [np.quantile(v[~np.isnan(v)], cover) if (~np.isnan(v)).sum() >= 10 else np.nan for v in scores.T]
    )
    c = pd.Series(c).interpolate(limit_direction="both").fillna(1.0)
    c = c.rolling(smooth, center=True, min_periods=1).median()
    return np.clip(c.to_numpy(), 0.25, 6.0)


#Evaluate the model by comparing the predicted values with the actual values. (CRPS)
def crps_ens(S, y):
    """S (n,H) ensemble, y (H,). Energy-form CRPS per horizon (NaN where y is NaN)."""
    n = S.shape[0]
    Ss = np.sort(S, axis=0)
    i = np.arange(1, n + 1)[:, None]
    return np.mean(np.abs(Ss - y[None, :]), axis=0) - np.sum((2 * i - n - 1) * Ss, axis=0) / n**2

#make the dataframe with the predicted values and the actual values. 
#It will also calculate the CRPS, coverage and width of the prediction interval.
def score_forecasts(Q, Y_raw, log, scale=None, n=199):
    rows = []
    for o in range(Q.shape[0]):
        s = quantiles_to_samples(Q[o], n)
        if scale is not None:
            s = apply_scale(s, scale)
        S = inv(s, log)
        y = Y_raw[o]
        lo, hi = np.quantile(S, [0.025, 0.975], axis=0)
        med = np.median(S, axis=0)
        ok = ~np.isnan(y)
        crps = crps_ens(S, np.where(ok, y, 0.0))
        for h in np.where(ok)[0]:
            rows.append((o, h + 1, y[h], med[h], crps[h], float(lo[h] <= y[h] <= hi[h]), hi[h] - lo[h]))
    d = pd.DataFrame(rows, columns=["origin", "h", "y", "med", "crps", "cover95", "width95"])
    d["abs_err"] = (d["y"] - d["med"]).abs()
    d["sq_err"] = (d["y"] - d["med"]) ** 2
    return d

#summarize the results by horizon or overall. 
#It will return a dataframe with the number of observations, MAE, RMSE, CRPS, coverage and width of the prediction interval.
def summarise(d, by_h=True):
    g = d.groupby("h") if by_h else d.assign(all=1).groupby("all")
    out = g.agg(n=("y", "size"), MAE=("abs_err", "mean"), RMSE=("sq_err", lambda v: np.sqrt(v.mean())),
                CRPS=("crps", "mean"), coverage95=("cover95", "mean"), width95=("width95", "mean"))
    out = out.reset_index(drop=not by_h)
    return out


#collect the context and the actual values for the origins.
def collect(series, backend, origins, ctx_len, log, horizon=HORIZON):
    ctxs, keep = [], []
    for o in origins:
        c = make_context(series, o, ctx_len, log)
        if c is not None:
            ctxs.append(c)
            keep.append(o)
    t0 = time.perf_counter()
    Q = backend.predict(ctxs, horizon)
    batch_s = time.perf_counter() - t0
    Y = np.full((len(keep), horizon), np.nan)
    for i, o in enumerate(keep):
        Y[i] = series.reindex(pd.date_range(o + pd.Timedelta(days=1), periods=horizon)).to_numpy()
    good = ~np.isnan(Y).all(axis=1)
    return Q[good], Y[good], batch_s / max(len(ctxs), 1), [k for k, g in zip(keep, good) if g]

#run all back tests
def cmd_evaluate(a):
    if not a.data:
       raise SystemExit("--data is required (path to the team's cleaned daily csv)")
    series = load_series(a.data, a.date_col, a.target)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    log = not a.no_log
    splits = {
        "val": pd.date_range(a.val_start, a.val_end, freq=f"{a.stride}D"),
        "test": pd.date_range(a.test_start, a.test_end, freq=f"{a.stride}D"),
    }
    by_h, overall, runtime = [], [], {}
    for name in a.models:
        t0 = time.perf_counter()
        backend = get_backend(name, a.ctx_len)
        load_s = time.perf_counter() - t0
        res = {}
        for sp, origins in splits.items():
            Q, Y, per_fc_s, kept = collect(series, backend, origins, a.ctx_len, log)
            res[sp] = (Q, Y)
            runtime[name] = {"load_s": round(load_s, 2), "per_35day_forecast_s": round(per_fc_s, 3)}
            print(f"[{name}] {sp}: {len(kept)} origins, {per_fc_s:.3f}s per 35-day forecast")
        scale = fit_scale(res["val"][0], fwd(res["val"][1], log), n=a.members)
        json.dump({"model": name, "log1p": log, "ctx_len": a.ctx_len, "scale": scale.tolist()},
                  open(out / f"calibration_{name}.json", "w"))
        for sp in splits:
            Q, Y = res[sp]
            variants = [("raw", None)] if sp == "val" else [("raw", None), ("calibrated", scale)]
            for tag, sc in variants:
                d = score_forecasts(Q, Y, log, sc, a.members)
                h = summarise(d).assign(model=name, split=sp, variant=tag)
                o = summarise(d, by_h=False).assign(model=name, split=sp, variant=tag)
                by_h.append(h)
                overall.append(o)
    by_h = pd.concat(by_h)
    overall = pd.concat(overall)
    by_h.to_csv(out / "metrics_by_horizon.csv", index=False)
    overall.to_csv(out / "metrics_overall.csv", index=False)
    json.dump(runtime, open(out / "runtime.json", "w"), indent=2)
    pd.set_option("display.width", 200)
    print("\n=== overall (all horizons pooled) ===")
    print(overall.drop(columns=["all"], errors="ignore").round(3).to_string(index=False))
    print("\n=== test, calibrated where available: selected horizons ===")
    sel = by_h[(by_h.split == "test") & (by_h.variant == "calibrated") & by_h.h.isin([1, 7, 14, 21, 28, 35])]
    print(sel[["model", "h", "n", "MAE", "RMSE", "CRPS", "coverage95", "width95"]].round(3).to_string(index=False))
    print(f"\nsaved to {out}/  (target coverage95 = 0.95, runtime budget = 20 min)")



# VERA forecast for submission
def to_vera(S, ref, model_id, site_id, depth_m, variable):
    n, H = S.shape
    dts = pd.date_range(ref + pd.Timedelta(days=1), periods=H)
    fmt = "%Y-%m-%d %H:%M:%S"
    return pd.DataFrame({
        "project_id": "vera4cast",
        "model_id": model_id,
        "datetime": np.tile(dts.strftime(fmt), n),
        "reference_datetime": ref.strftime(fmt),
        "duration": "P1D",
        "site_id": site_id,
        "depth_m": depth_m,
        "family": "ensemble",
        "parameter": np.repeat(np.arange(1, n + 1), H),
        "variable": variable,
        "prediction": S.ravel(),
    })

# actual forecast for submission (CMD)
def cmd_forecast(a):
    if not a.data:
       raise SystemExit("--data is required (path to the team's cleaned daily csv)")
    series = load_series(a.data, a.date_col, a.target)
    log = not a.no_log
    ref = pd.Timestamp(a.reference_date) if a.reference_date else pd.Timestamp.utcnow().tz_localize(None).normalize()
    last = series.loc[:ref].last_valid_index()
    lag = (ref - last).days
    ctx = make_context(series, last, a.ctx_len, log, max_stale=a.max_stale)
    if ctx is None:
        raise SystemExit(f"No usable context: last observation {last.date()} is {lag} days before {ref.date()}")
    t0 = time.perf_counter()
    backend = get_backend(a.model, a.ctx_len)
    load_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    total = HORIZON + lag
    samples = quantiles_to_samples(backend.predict([ctx], total)[0], a.members)
    if a.calib:
        c = np.array(json.load(open(a.calib))["scale"])
        c = np.r_[c, np.full(max(0, total - len(c)), c[-1])][:total]
        samples = apply_scale(samples, c)
    S = inv(samples, log)[:, lag: lag + HORIZON]      # keep ref+1 ... ref+35
    df = to_vera(S, ref, a.model_id or f"chlorocast_{a.model}", a.site, a.depth, a.target)
    assert df["prediction"].notna().all() and (df["prediction"] >= 0).all()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    f = out / f"daily-{ref.date()}-{df['model_id'].iat[0]}.csv.gz"
    df.to_csv(f, index=False)
    print(f"wrote {f} ({len(df)} rows); context ends {last.date()} (lag {lag}d); "
          f"load {load_s:.1f}s, forecast {time.perf_counter() - t0:.2f}s")


# run the function
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(s):
        s.add_argument("--data", help="cleaned daily csv; omit to pull targets from VERA")
        s.add_argument("--date-col")
        s.add_argument("--target", default=VARIABLE)
        s.add_argument("--site", default="fcre")
        s.add_argument("--depth", type=float, default=1.6)
        s.add_argument("--ctx-len", type=int, default=512)
        s.add_argument("--members", type=int, default=199, help="ensemble size (odd)")
        s.add_argument("--no-log", action="store_true", help="skip the log1p transform")

    e = sub.add_parser("evaluate")
    common(e)
    e.add_argument("--models", nargs="+", default=["persistence", "chronos", "timesfm"])
    e.add_argument("--stride", type=int, default=3, help="days between forecast origins")
    e.add_argument("--val-start", default="2025-01-01")
    e.add_argument("--val-end", default="2025-12-31")
    e.add_argument("--test-start", default="2026-01-01")
    e.add_argument("--test-end", default="2026-09-20")
    e.add_argument("--out", default="results")
    e.set_defaults(fn=cmd_evaluate)

    f = sub.add_parser("forecast")
    common(f)
    f.add_argument("--model", required=True, choices=["persistence", "chronos", "timesfm"])
    f.add_argument("--model-id")
    f.add_argument("--calib", help="results/calibration_<model>.json from evaluate")
    f.add_argument("--reference-date", help="YYYY-MM-DD, default today (UTC)")
    f.add_argument("--max-stale", type=int, default=7)
    f.add_argument("--out", default="forecasts")
    f.set_defaults(fn=cmd_forecast)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()