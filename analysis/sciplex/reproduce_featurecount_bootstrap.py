#!/usr/bin/env python3
"""Reproduce sci-Plex 500/1000/2000-gene robustness and drug-level bootstrap.

The input matrix must contain the nested top-2000 variance-ranked genes in order.
Bootstrap contract: 20,000 samples per feature count; 188 drugs sampled with
replacement per sample; one NumPy default_rng(20260913) used sequentially for
500, 1000, then 2000 genes. CIs are empirical 2.5/97.5 percentiles and selection
frequency is the fraction of samples in which a resolution has minimum SSE.
"""
import argparse
import numpy as np
import pandas as pd

LEVELS=["Shared","Cell line","Cell line × 2 dose bins","Cell line × exact dose"]
NGENES=[500,1000,2000]; SEED=20260913; B=20000

def context_cols(level):
    return {"Shared":[],"Cell line":["cell_line"],"Cell line × 2 dose bins":["cell_line","dose_bin"],"Cell line × exact dose":["cell_line","dose"]}[level]

def wm(X,w): return np.average(X,axis=0,weights=np.asarray(w,float))

def per_drug_sse(meta,R,level,drugs):
    cols=context_cols(level); out=np.zeros(len(drugs),float); reps=sorted(meta.replicate.astype(str).unique())
    for j,drug in enumerate(drugs):
        idx=np.flatnonzero(meta.perturbation.astype(str).to_numpy()==drug); total=0.0
        for held in reps:
            te=idx[meta.iloc[idx].replicate.astype(str).to_numpy()==held]
            tr=idx[meta.iloc[idx].replicate.astype(str).to_numpy()!=held]
            if not cols: pred={():wm(R[tr],meta.iloc[tr].n_match)}; keys=[()]*len(te)
            else:
                pred={}
                for key,g in meta.iloc[tr].groupby(cols,sort=True):
                    key=key if isinstance(key,tuple) else (key,); ii=g.index.to_numpy(); pred[key]=wm(R[ii],g.n_match)
                keys=[tuple(row[c] for c in cols) for _,row in meta.iloc[te].iterrows()]
            for ii,key in zip(te,keys): d=R[ii]-pred[key]; total+=float(d@d)
        out[j]=total
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--meta',required=True); ap.add_argument('--matrix',required=True)
    ap.add_argument('--summary-out',default='SciPlex_featurecount_robustness_summary.csv')
    ap.add_argument('--bootstrap-out',default='SciPlex_featurecount_robustness_bootstrap.csv'); args=ap.parse_args()
    meta=pd.read_csv(args.meta).copy(); Rall=np.asarray(np.load(args.matrix)['response'],float)
    meta['dose_bin']=np.where(meta.dose.astype(float)<=100.0,'low','high')
    drugs=np.array(sorted(meta.perturbation.astype(str).unique())); n=len(drugs)
    if n!=188: raise ValueError(f'Expected 188 drugs, found {n}')
    rng=np.random.default_rng(SEED); summary=[]; boots=[]
    for ng in NGENES:
        R=Rall[:,:ng]; pdse=np.vstack([per_drug_sse(meta,R,l,drugs) for l in LEVELS])
        totals=pdse.sum(axis=1); base=totals[0]
        for k,l in enumerate(LEVELS): summary.append((ng,l,1-totals[k]/base,totals[k],int(np.argsort(totals).tolist().index(k)+1)))
        gains=np.empty((B,len(LEVELS))); best=np.zeros((B,len(LEVELS)),bool)
        for b in range(B):
            idx=rng.choice(n,n,replace=True); ss=pdse[:,idx].sum(axis=1); gains[b]=1-ss/ss[0]; best[b,np.argmin(ss)]=True
        for k,l in enumerate(LEVELS):
            boots.append((ng,l,float(gains[:,k].mean()),float(np.percentile(gains[:,k],2.5)),float(np.percentile(gains[:,k],97.5)),float(best[:,k].mean())))
    pd.DataFrame(summary,columns=['n_genes','resolution','gain_vs_shared','total_SSE','rank']).to_csv(args.summary_out,index=False)
    out=pd.DataFrame(boots,columns=['n_genes','resolution','boot_mean_gain','CI_low','CI_high','P_best'])
    out.to_csv(args.bootstrap_out,index=False); print(out.to_string(index=False))
if __name__=='__main__': main()
