from __future__ import annotations
import numpy as np
import pandas as pd
from numba import njit


@njit(cache=True)
def _group_sufficient(
    X,
    src,
    parent,
    child,
    weight,
    n_sources,
    n_parents,
    n_children,
):
    n,p=X.shape

    # Weighted sums by parent / source-parent / child / source-child.
    sp=np.zeros((n_parents,p))
    wsp=np.zeros(n_parents)
    ssp=np.zeros((n_sources,n_parents,p))
    wssp=np.zeros((n_sources,n_parents))

    sc=np.zeros((n_children,p))
    wsc=np.zeros(n_children)
    ssc=np.zeros((n_sources,n_children,p))
    wssc=np.zeros((n_sources,n_children))

    for i in range(n):
        s=src[i]
        pa=parent[i]
        ch=child[i]
        w=weight[i]
        wsp[pa]+=w
        wssp[s,pa]+=w
        wsc[ch]+=w
        wssc[s,ch]+=w
        for f in range(p):
            v=w*X[i,f]
            sp[pa,f]+=v
            ssp[s,pa,f]+=v
            sc[ch,f]+=v
            ssc[s,ch,f]+=v

    A=np.zeros(n_sources)
    B=np.zeros(n_sources)
    C=np.zeros(n_sources)

    for i in range(n):
        s=src[i]
        pa=parent[i]
        ch=child[i]

        dc=wsp[pa]-wssp[s,pa]
        df=wsc[ch]-wssc[s,ch]
        if dc<=0.0 or df<=0.0:
            continue

        Ai=0.0
        Bi=0.0
        Ci=0.0
        for f in range(p):
            muc=(sp[pa,f]-ssp[s,pa,f])/dc
            muf=(sc[ch,f]-ssc[s,ch,f])/df
            a=X[i,f]-muc
            c=muf-muc
            Ai += a*a
            Bi += a*c
            Ci += c*c

        A[s]+=Ai
        B[s]+=Bi
        C[s]+=Ci

    return A,B,C


def prepare_fast_common_taxonomy(response, meta, parent_col, child_col, weight_col):
    X=np.asarray(response,dtype=np.float64)
    m=meta.reset_index(drop=True).copy()

    units=pd.Index(pd.unique(m["biological_unit"]))
    src_map={u:i for i,u in enumerate(units)}
    src=m["biological_unit"].map(src_map).to_numpy(np.int64)

    # Parent code.
    parents=pd.Index(pd.unique(m[parent_col]))
    parent_map={u:i for i,u in enumerate(parents)}
    parent=m[parent_col].map(parent_map).to_numpy(np.int64)

    # Child labels are parent-qualified so identical strings under different
    # parents can never be mixed.
    child_key=m[parent_col].astype(str)+"<<<>>>"+m[child_col].astype(str)
    child_values=pd.Index(pd.unique(child_key))
    child_map={u:i for i,u in enumerate(child_values)}
    child=child_key.map(child_map).to_numpy(np.int64)

    if weight_col is None:
        weight=np.ones(len(m),dtype=np.float64)
    else:
        weight=m[weight_col].to_numpy(np.float64)

    # Domains for common relabeling: one child-label set per source × parent.
    domains={}
    for s in range(len(units)):
        for pa in range(len(parents)):
            ix=np.where((src==s)&(parent==pa))[0]
            if len(ix):
                domains[(s,pa)]=np.unique(child[ix])

    # Split by perturbation. Cross-perturbation inner products are never used.
    groups=[]
    for pert,g in m.groupby("perturbation",sort=False):
        ix=g.index.to_numpy(np.int64)
        groups.append(ix)

    return {
        "X":X,
        "meta":m,
        "src":src,
        "parent":parent,
        "child":child,
        "weight":weight,
        "domains":domains,
        "groups":groups,
        "n_sources":len(units),
        "n_parents":len(parents),
        "n_children":len(child_values),
    }


def permute_child_codes(prepared, rng):
    child=prepared["child"]
    src=prepared["src"]
    parent=prepared["parent"]
    mapped=child.copy()

    for (s,pa),vals in prepared["domains"].items():
        perm=rng.permutation(vals)
        # mapping shared across all perturbations
        for a,b in zip(vals,perm):
            ix=(src==s)&(parent==pa)&(child==a)
            mapped[ix]=b
    return mapped


def sufficient_from_codes(prepared, mapped_child):
    A=np.zeros(prepared["n_sources"])
    B=np.zeros(prepared["n_sources"])
    C=np.zeros(prepared["n_sources"])

    for ix in prepared["groups"]:
        a,b,c=_group_sufficient(
            prepared["X"][ix],
            prepared["src"][ix],
            prepared["parent"][ix],
            mapped_child[ix],
            prepared["weight"][ix],
            prepared["n_sources"],
            prepared["n_parents"],
            prepared["n_children"],
        )
        A+=a; B+=b; C+=c

    return A,B,C


def nested_partial_gain(A,B,C):
    At=A.sum()
    Bt=B.sum()
    Ct=C.sum()

    ws=np.zeros(len(A))
    gains=np.zeros(len(A))

    for s in range(len(A)):
        trainB=Bt-B[s]
        trainC=Ct-C[s]
        w=0.0 if trainC<=0 else float(np.clip(trainB/trainC,0.0,1.0))
        ws[s]=w
        gains[s]=2*w*B[s]-w*w*C[s]

    raw=np.nan if At<=0 else float(np.sum(2*B-C)/At)
    partial=np.nan if At<=0 else float(gains.sum()/At)

    return raw,partial,float(ws.mean())
