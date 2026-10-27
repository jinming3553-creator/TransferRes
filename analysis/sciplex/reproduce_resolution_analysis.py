#!/usr/bin/env python3
"""Reproduce the primary sci-Plex resolution-selection analysis.

Inputs are the response matrix and metadata produced by build_response_matrix.py.
For each candidate resolution, prediction is evaluated by reciprocal replicate
holdout. Training-context centroids are weighted by n_match and raw squared
Euclidean errors are summed over held-out units and genes.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

LEVELS = ["Shared", "Cell line", "Cell line × 2 dose bins", "Cell line × exact dose"]

def context_cols(level):
    return {"Shared": [], "Cell line": ["cell_line"],
            "Cell line × 2 dose bins": ["cell_line", "dose_bin"],
            "Cell line × exact dose": ["cell_line", "dose"]}[level]

def weighted_mean(X, w):
    return np.average(X, axis=0, weights=np.asarray(w, float))

def total_sse(meta, R, level):
    cols=context_cols(level); total=0.0
    reps=sorted(meta.replicate.astype(str).unique())
    if len(reps)!=2: raise ValueError(f"Expected 2 replicates, found {reps}")
    for drug in sorted(meta.perturbation.astype(str).unique()):
        idx=np.flatnonzero(meta.perturbation.astype(str).to_numpy()==drug)
        for held in reps:
            te=idx[meta.iloc[idx].replicate.astype(str).to_numpy()==held]
            tr=idx[meta.iloc[idx].replicate.astype(str).to_numpy()!=held]
            if not cols:
                pred={(): weighted_mean(R[tr], meta.iloc[tr].n_match)}; keys=[()]*len(te)
            else:
                pred={}
                for key,g in meta.iloc[tr].groupby(cols, sort=True):
                    key=key if isinstance(key,tuple) else (key,)
                    ii=g.index.to_numpy(); pred[key]=weighted_mean(R[ii],g.n_match)
                keys=[tuple(row[c] for c in cols) for _,row in meta.iloc[te].iterrows()]
            for ii,key in zip(te,keys):
                d=R[ii]-pred[key]; total += float(d@d)
    return total

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--meta', required=True)
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--n-genes', type=int, default=1000)
    ap.add_argument('--out', default='SciPlex_resolution_summary.csv')
    args=ap.parse_args()
    meta=pd.read_csv(args.meta).copy(); z=np.load(args.matrix); R=np.asarray(z['response'],float)
    if 'unit_id' in z and not np.array_equal(np.asarray(z['unit_id']),meta.unit_id.to_numpy()):
        raise ValueError('Matrix and metadata unit_id are not aligned')
    if args.n_genes>R.shape[1]: raise ValueError('Requested more genes than matrix contains')
    R=R[:,:args.n_genes]
    meta['dose_bin']=np.where(meta.dose.astype(float)<=100.0,'low','high')
    ss={level:total_sse(meta,R,level) for level in LEVELS}; base=ss['Shared']
    out=pd.DataFrame({'resolution':LEVELS,'total_SSE':[ss[x] for x in LEVELS]})
    out['gain_vs_shared']=1-out.total_SSE/base
    out['rank']=out.total_SSE.rank(method='min').astype(int)
    out.to_csv(args.out,index=False)
    print(out.to_string(index=False)); print('Selected:',out.loc[out.total_SSE.idxmin(),'resolution'])
if __name__=='__main__': main()
