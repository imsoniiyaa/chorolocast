"""Causal, anchored analog distributions for ChloroFM experiments and forecasting."""
from __future__ import annotations

import numpy as np
import pandas as pd
from chlorofm import TARGET, make_features


def analog_bank(d, origins, members=199):
    """Build marginal ensembles from fully observed historical trajectories.

    Missing future labels are never imputed. Historical features are forward
    filled; scaling and neighbor selection use only history available at origin.
    Anchoring transfers historical changes to the current chlorophyll level,
    with a decay back to the historical level for longer forecast horizons.
    """
    f = make_features(d).replace([np.inf, -np.inf], np.nan)
    # No future backfill, including at the beginning of the record.
    groups = {
        'season': ['sin', 'cos'],
        'state': ['sin', 'cos', 'ch_last', 'ch_mean7', 'ch_mean30', 'ch_slope7'],
        'environment': list(f.columns),
    }
    raw = d[TARGET].to_numpy(float)
    paths = np.lib.stride_tricks.sliding_window_view(raw, 36)
    valid = np.isfinite(paths).all(axis=1)
    u = (np.arange(members) + .5) / members
    bank = {}
    for oi, origin in enumerate(origins):
        pos = d.index.get_loc(origin)
        eligible = np.flatnonzero(valid & (np.arange(len(paths)) >= 60)
                                  & (np.arange(len(paths)) + 35 <= pos))
        if len(eligible) < 120:
            raise ValueError('need at least 120 complete historical trajectories')
        current = float(d[TARGET].iloc[:pos + 1].ffill().iloc[-1])
        for group, cols in groups.items():
            x = f[cols].iloc[eligible].to_numpy(float)
            scale = np.nanstd(x, axis=0, ddof=1)
            scale = np.where(np.isfinite(scale) & (scale > 1e-8), scale, 1.)
            delta = (x - f[cols].iloc[pos].to_numpy(float)) / scale
            dist = np.nanmean(delta ** 2, axis=1)
            ranking = eligible[np.argsort(np.nan_to_num(dist, nan=np.inf), kind='stable')]
            for k in (30, 60, 120):
                trajectories = paths[ranking[:k]]
                for tau in (0, 7, 21, 60):
                    decay = np.zeros(35) if tau == 0 else np.exp(-np.arange(1, 36) / tau)
                    for mode in ('raw', 'log') if tau else ('raw',):
                        future = trajectories[:, 1:]
                        if mode == 'log':
                            adjusted = np.expm1(np.log1p(future) + decay *
                                (np.log1p(current) - np.log1p(trajectories[:, :1])))
                        else:
                            adjusted = future + decay * (current - trajectories[:, :1])
                        name = f'{group}_k{k}_{mode}_tau{tau}'
                        if name not in bank:
                            bank[name] = np.empty((len(origins), members, 35), np.float32)
                        bank[name][oi] = np.maximum(np.quantile(adjusted, u, axis=0), 0)
    return bank


def distribution_mix(a, b, weight):
    """Approximate a mixture CDF on equal-mass quantiles, retaining disagreement."""
    if weight <= 0:
        return b.copy()
    if weight >= 1:
        return a.copy()
    out = np.empty_like(a)
    n = a.shape[1]
    u = (np.arange(n) + .5) / n
    for i in range(len(a)):
        for h in range(a.shape[2]):
            values = np.r_[a[i, :, h], b[i, :, h]]
            mass = np.r_[np.full(n, weight / n), np.full(n, (1 - weight) / n)]
            order = np.argsort(values, kind='stable')
            out[i, :, h] = np.interp(u, np.cumsum(mass[order]) - mass[order] / 2, values[order])
    return out
