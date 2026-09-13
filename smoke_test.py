
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

from transferres import (
    fit_resolution,
    bootstrap_biological_units,
    diagnose_identifiability,
)

meta=pd.read_csv(ROOT/"sample_data"/"sample_metadata.csv")
X=np.load(ROOT/"sample_data"/"sample_response.npy")
meta=meta.rename(columns={"source":"biological_unit","gene_target":"perturbation"})

fit=fit_resolution(
    X,meta,
    hierarchy=["broad","class_name"],
    weight_col="n_match"
)
assert fit.selected_level=="class_name"

boot=bootstrap_biological_units(fit,n_boot=50,seed=1)
ident=diagnose_identifiability(meta,parent_col="broad",child_col="class_name")

print("PASS",fit.selected_level,int(ident.summary_table.iloc[0].n_biological_units_total))
