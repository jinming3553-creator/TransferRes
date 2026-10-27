#!/usr/bin/env python3
"""Build the frozen sci-Plex response representation from the source h5ad.

Population: A549, K562 and MCF7 cells at 24 h. Genes are ranked once by sample
variance across this population; the first 2000 form the robustness matrix and
the first 1000 are the primary feature set. Perturbation responses are group
means minus vehicle-control means matched within replicate x cell line. Vehicle
means are cell-count weighted and n_match=min(n_pert,n_ctrl).
"""
import argparse, hashlib
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--h5ad',required=True); ap.add_argument('--outdir',default='.')
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    a=ad.read_h5ad(args.h5ad,backed='r'); obs=a.obs.copy().reset_index(drop=True)
    mask=obs.cell_line.astype(str).isin(['A549','K562','MCF7']) & (obs.time.astype(float)==24.0)
    idx=np.flatnonzero(mask.to_numpy()); p=a.n_vars; s=np.zeros(p); ss=np.zeros(p); n=0
    for z in range(0,len(idx),2048):
        q=idx[z:z+2048]; X=a.X[q,:]
        if sparse.issparse(X): s+=np.asarray(X.sum(0)).ravel(); ss+=np.asarray(X.multiply(X).sum(0)).ravel()
        else: X=np.asarray(X,float); s+=X.sum(0); ss+=(X*X).sum(0)
        n+=len(q)
    var=(ss-s*s/n)/max(n-1,1); rank=np.argsort(-var,kind='stable')[:2000]
    features=pd.DataFrame({'rank':np.arange(1,2001),'feature_index':rank,'gene':a.var_names[rank].astype(str),'variance':var[rank]})
    features.to_csv(out/'SciPlex_features_top2000_nested.csv',index=False)
    o=obs.iloc[idx].copy().reset_index(drop=True)
    for c in ['replicate','cell_line','perturbation','dose_value']: o[c]=o[c].astype(str)
    gcols=['replicate','cell_line','dose_value','perturbation']; mi=pd.MultiIndex.from_frame(o[gcols].astype('string').fillna('<NA>'))
    codes,uniques=pd.factorize(mi,sort=True); sums=np.zeros((len(uniques),2000)); counts=np.zeros(len(uniques),np.int64)
    for z in range(0,len(idx),4096):
        q=idx[z:z+4096]; X=a.X[q,:][:,rank]
        if sparse.issparse(X): X=X.toarray()
        X=np.asarray(X,float); cc=codes[z:z+len(q)]; np.add.at(sums,cc,X); np.add.at(counts,cc,1)
    groups=uniques.to_frame(index=False); groups.columns=gcols; groups['n_cells']=counts; means=sums/counts[:,None]
    ctrl=groups.perturbation.eq('control'); ctrlmap={}
    for (rep,cl),g in groups[ctrl].groupby(['replicate','cell_line'],observed=True):
        ii=g.index.to_numpy(); w=groups.loc[ii,'n_cells'].to_numpy(float); ctrlmap[(str(rep),str(cl))]=(np.average(means[ii],axis=0,weights=w),int(w.sum()))
    rows=[]; vec=[]
    for i,r in groups[~ctrl].iterrows():
        key=(str(r.replicate),str(r.cell_line))
        if key not in ctrlmap: continue
        cm,nc=ctrlmap[key]; rows.append({'perturbation':str(r.perturbation),'replicate':str(r.replicate),'cell_line':str(r.cell_line),'dose':float(r.dose_value),'time':24.0,'n_pert':int(r.n_cells),'n_ctrl':nc,'n_match':min(int(r.n_cells),nc)}); vec.append((means[i]-cm).astype(np.float32))
    meta=pd.DataFrame(rows).sort_values(['perturbation','replicate','cell_line','dose']).reset_index(drop=True); meta.insert(0,'unit_id',np.arange(1,len(meta)+1)); R=np.vstack(vec)
    # Reorder R consistently with sorted metadata.
    # Reconstruct row order mapping from the unsorted rows to avoid dependence on group iteration order.
    uns=pd.DataFrame(rows); uns['_row']=np.arange(len(uns)); keys=['perturbation','replicate','cell_line','dose','time']
    order=meta[keys].merge(uns[keys+['_row']],on=keys,validate='one_to_one')['_row'].to_numpy(int); R=R[order]
    meta.to_csv(out/'SciPlex_featurecount_meta.csv',index=False)
    np.savez_compressed(out/'SciPlex_response_matrix_top2000_nested.npz',response=R,feature_index=rank.astype(np.int64))
    R1=R[:,:1000]; np.savez_compressed(out/'SciPlex_response_matrix_1000genes.npz',response=R1,unit_id=meta.unit_id.to_numpy(np.int64),feature_index=rank[:1000].astype(np.int64))
    features.iloc[:1000].to_csv(out/'SciPlex_response_matrix_features.csv',index=False)
    sha=hashlib.sha256((out/'SciPlex_response_matrix_1000genes.npz').read_bytes()).hexdigest()
    pd.DataFrame([{'n_units':len(meta),'n_features':1000,'n_drugs':meta.perturbation.nunique(),'n_replicates':meta.replicate.nunique(),'n_cell_lines':meta.cell_line.nunique(),'n_doses':meta.dose.nunique(),'time':24.0,'npz_sha256':sha}]).to_csv(out/'SciPlex_response_matrix_QC.csv',index=False)
    a.file.close(); print(f'Built {len(meta)} units; {meta.perturbation.nunique()} drugs; response shape {R1.shape}')
if __name__=='__main__': main()
