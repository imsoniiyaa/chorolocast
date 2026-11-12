"""Experimental equal blend of ChloroFM v1 and seasonal-residual TimesFM.

This design was informed by 2026 development backtests. Evaluation on those
dates is descriptive, not an independent test of generalization.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from chlorofm import (TARGET, load_daily, load_config, foundation_samples,
                      hybrid_samples, outcomes, crps_ensemble, TimesFM25, to_vera)
from chlorofm_seasonal import seasonal_foundation


def predict(d, dates, config, members=199, baseline=False):
    backend=TimesFM25(1024)
    foundation=foundation_samples(backend,d[TARGET],dates,members,1024)
    original=hybrid_samples(d,dates,foundation,config,members)
    seasonal=seasonal_foundation(backend,d[TARGET],dates,1024,.5,members)
    result={'chlorofm_consensus':.5*original+.5*seasonal,
            'chlorofm_v1':original,'timesfm_log1024':foundation}
    if baseline:
        result['timesfm_log512']=foundation_samples(backend,d[TARGET],dates,members,512)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['evaluate','forecast'])
    p.add_argument('--data',required=True)
    p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config/chlorofm.json'))
    p.add_argument('--start',default='2026-01-01')
    p.add_argument('--end',default='2026-08-16')
    p.add_argument('--stride',type=int,default=7)
    p.add_argument('--reference-date')
    p.add_argument('--out',default='results_chlorofm_consensus')
    p.add_argument('--site',default='fcre')
    p.add_argument('--depth',type=float,default=1.6)
    a=p.parse_args()
    d=load_daily(a.data); config=load_config(a.config)
    out=Path(a.out); out.mkdir(exist_ok=True,parents=True)
    if a.command=='forecast':
        ref=pd.Timestamp(a.reference_date) if a.reference_date else d[TARGET].last_valid_index()
        if ref not in d.index or not np.isfinite(d.loc[ref,TARGET]):
            raise ValueError('forecast reference must have an observed target; update the input first')
        results=predict(d,pd.DatetimeIndex([ref]),config)
        frame=to_vera(results['chlorofm_consensus'][0],ref,'chlorocast_consensus_experimental',a.site,a.depth,TARGET)
        path=out/f'daily-{ref.date()}-chlorocast_consensus_experimental.csv.gz'
        frame.to_csv(path,index=False); print(path)
        return
    if a.stride<1: raise ValueError('stride must be positive')
    dates=pd.date_range(a.start,a.end,freq=f'{a.stride}D')
    if len(dates)==0 or dates[-1]+pd.Timedelta(days=35)>d.index[-1]:
        raise ValueError('evaluation must end at least 35 days before data ends')
    results=predict(d,dates,config,baseline=True)
    y=outcomes(d,dates); rows=[]
    for name,samples in results.items():
        c=crps_ensemble(samples,y)
        lower,upper=np.quantile(samples,[.025,.975],axis=1)
        for i,origin in enumerate(dates):
            for h in range(35):
                if np.isfinite(y[i,h]):
                    rows.append(dict(model=name,origin=str(origin.date()),h=h+1,
                        CRPS=c[i,h],coverage95=float(lower[i,h]<=y[i,h]<=upper[i,h]),
                        width95=upper[i,h]-lower[i,h]))
    detail=pd.DataFrame(rows)
    detail.to_csv(out/'scores.csv',index=False)
    overall=detail.groupby('model').agg(CRPS=('CRPS','mean'),n=('CRPS','size'),
                                        coverage95=('coverage95','mean'),width95=('width95','mean'))
    overall.to_csv(out/'metrics_overall.csv')
    detail.groupby(['model','h']).agg(CRPS=('CRPS','mean'),n=('CRPS','size')).to_csv(out/'metrics_by_horizon.csv')
    meta=dict(evaluation_type='development backtest; 2026 informed architecture selection',
              data_sha256=hashlib.sha256(Path(a.data).read_bytes()).hexdigest(),
              start=a.start,end=a.end,stride=a.stride,origins=len(dates),members=199,
              model='50% ChloroFM v1 + 50% TimesFM1024 with half seasonal removal')
    (out/'provenance.json').write_text(json.dumps(meta,indent=2))
    print(overall.round(5).to_string())


if __name__=='__main__': main()
