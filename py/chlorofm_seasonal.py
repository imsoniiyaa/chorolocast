"""Foundation forecasts of deviations from a causally estimated annual cycle."""
import numpy as np
import pandas as pd
from fm import quantiles_to_samples


def design(dates):
    day=dates.dayofyear.to_numpy()/365.25
    return np.stack([np.ones(len(day)), *[fun(2*np.pi*k*day)
                    for k in (1,2,3) for fun in (np.sin,np.cos)]],axis=1)


def seasonal_foundation(backend,series,origins,context=1024,strength=1.,members=199):
    contexts=[]; seasonal=[]
    for origin in origins:
        history=series.loc[:origin].iloc[-1461:]
        valid=history.notna()
        X=design(history.index[valid]); y=np.log1p(np.maximum(history[valid].to_numpy(float),0))
        penalty=np.eye(X.shape[1])*10; penalty[0,0]=1e-6
        coef=np.linalg.solve(X.T@X+penalty,X.T@y)
        ctx=history.iloc[-context:].interpolate(limit_direction='both')
        contexts.append((np.log1p(np.maximum(ctx.to_numpy(float),0))-strength*design(ctx.index)@coef).astype(np.float32))
        future=pd.date_range(origin+pd.Timedelta(days=1),periods=35)
        seasonal.append(strength*design(future)@coef)
    Q=backend.predict(contexts,35)
    return np.stack([np.maximum(np.expm1(np.minimum(quantiles_to_samples(q,members)+season,15)),0)
                     for q,season in zip(Q,seasonal)])
