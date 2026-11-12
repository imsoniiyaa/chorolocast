"""ens4: the exact fixed ChloroFM seasonal consensus (199 raw-scale members)."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from weakref import ref
import numpy as np
import pandas as pd
from chlorofm import TARGET, load_daily, load_config, outcomes, crps_ensemble, to_vera
from chlorofm_consensus import predict


def score(samples, y, dates):
    crps = crps_ensemble(samples, y)
    low, high = np.quantile(samples, [.025, .975], axis=1)
    med = np.median(samples, axis=1)
    rows = []
    for i, date in enumerate(dates):
        for h in range(35):
            if np.isfinite(y[i, h]):
                error = med[i, h] - y[i, h]
                rows.append(dict(origin=str(date.date()), h=h+1, y=y[i,h], med=med[i,h],
                    crps=crps[i,h], cover95=float(low[i,h] <= y[i,h] <= high[i,h]),
                    width95=high[i,h]-low[i,h], abs_err=abs(error), sq_err=error**2))
    return pd.DataFrame(rows)


def summary(d):
    return dict(n=len(d), MAE=d.abs_err.mean(), RMSE=np.sqrt(d.sq_err.mean()),
        CRPS=d.crps.mean(), coverage95=d.cover95.mean(), width95=d.width95.mean())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['forecast', 'evaluate'])
    p.add_argument('--data', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--config', default=str(Path(__file__).resolve().parents[1]/'config/chlorofm.json'))
    p.add_argument('--reference-date')
    p.add_argument('--model-id', default='chlorocast_ens4')
    p.add_argument('--site', default='fcre')
    p.add_argument('--depth', type=float, default=1.6)
    p.add_argument('--val-start', default='2025-01-01')
    p.add_argument('--val-end', default='2025-12-31')
    p.add_argument('--test-start', default='2026-01-01')
    p.add_argument('--test-end', default='2026-09-20')
    p.add_argument('--stride', type=int, default=3)
    a = p.parse_args()
    d = load_daily(a.data); config = load_config(a.config)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.command == 'forecast':
        ref = pd.Timestamp(a.reference_date) if a.reference_date else d[TARGET].last_valid_index()
        if ref not in d.index:
            raise ValueError('ens4 requires an observed target on the reference date')
        if not np.isfinite(d[TARGET].ffill().loc[ref]):
            raise ValueError('ens4 requires an observed target on the reference date')
        samples = predict(d, pd.DatetimeIndex([ref]), config)['chlorofm_consensus'][0]
        frame = to_vera(samples, ref, a.model_id, a.site, a.depth, TARGET)
        assert len(frame) == 6965 and np.isfinite(frame.prediction).all() and (frame.prediction >= 0).all()
        path = out/f'{a.model_id}.csv'; frame.to_csv(path, index=False)
        print(f'wrote {path} ({len(frame)} rows)')
        return
    if a.stride < 1:
        raise ValueError('stride must be positive')
    overall, by_h, periods = [], [], {}
    start_time = time.perf_counter(); count = 0
    for split in ['val', 'test']:
        start = pd.Timestamp(getattr(a, split+'_start'))
        end = min(pd.Timestamp(getattr(a, split+'_end')), d.index[-1])
        # Keep every scored target inside its split; no end-of-validation spillover.
        dates = pd.date_range(start, end-pd.Timedelta(days=35), freq=f'{a.stride}D')
        if not len(dates):
            raise ValueError(f'{split} needs at least 36 days')
        samples = predict(d, dates, config)['chlorofm_consensus']
        detail = score(samples, outcomes(d, dates), dates)
        detail.to_csv(out/f'scores_ens4_{split}_raw.csv', index=False)
        tags = dict(model='ens4', split=split, variant='raw')
        overall.append(dict(**tags, **summary(detail)))
        for h, rows in detail.groupby('h'):
            by_h.append(dict(**tags, h=h, **summary(rows)))
        periods[split] = dict(first_origin=str(dates[0].date()), last_origin=str(dates[-1].date()), origins=len(dates))
        count += len(dates)
    pd.DataFrame(overall).to_csv(out/'metrics_overall.csv', index=False)
    pd.DataFrame(by_h).to_csv(out/'metrics_by_horizon.csv', index=False)
    pd.DataFrame([dict(model='ens4', load_s=0, per_35day_forecast_s=(time.perf_counter()-start_time)/count)]).to_csv(out/'runtime.csv', index=False)
    (out/'provenance.json').write_text(json.dumps(dict(model='ens4', members=199,
        data_sha256=hashlib.sha256(Path(a.data).read_bytes()).hexdigest(), periods=periods,
        evaluation_type='development backtest; prior architecture selection used 2026',
        calibration='none; exact fixed raw-sample consensus'), indent=2))
    print(pd.DataFrame(overall).to_string(index=False))


if __name__ == '__main__':
    main()
