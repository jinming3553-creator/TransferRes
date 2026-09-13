
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

from transferres import (
    fit_resolution,
    bootstrap_biological_units,
    diagnose_identifiability,
    diagnose_metric_sensitivity,
    diagnose_feature_count_sensitivity,
    common_taxonomy_permutation_test,
    calibrate_signal_detectability,
)
from interpretation import build_interpretation

meta=pd.read_csv(ROOT/"sample_data"/"sample_metadata.csv")
X=np.load(ROOT/"sample_data"/"sample_response.npy")
meta=meta.rename(columns={"source":"biological_unit","gene_target":"perturbation"})
hier=["broad","class_name"]

report={}

# CORE
fit=fit_resolution(X,meta,hierarchy=hier,weight_col="n_match")
assert fit.selected_level=="class_name"
report["core"]="PASS"

# RECOMMENDED
boot=bootstrap_biological_units(fit,n_boot=100,seed=1)
ident=diagnose_identifiability(meta,parent_col="broad",child_col="class_name")
assert not boot.selection_frequency.empty
assert not ident.support_frontier.empty
report["recommended"]="PASS"

# FULL
metric=diagnose_metric_sensitivity(
    X,meta,hier,weight_col="n_match",
    primary_level=fit.selected_level,n_boot=100,seed=1
)

# requested GUI defaults 500/1000/2000 must safely collapse to available dimension
requested=[500,1000,2000]
p=X.shape[1]
counts=sorted({min(x,p) for x in requested})
fs=diagnose_feature_count_sensitivity(
    X,meta,hier,
    feature_counts=counts,
    weight_col="n_match",
    reference_level=fit.selected_level
)

perm=common_taxonomy_permutation_test(
    X,meta,parent_col="broad",child_col="class_name",
    weight_col="n_match",n_permutations=49,seed=4
)

ft=fit.fold_table
q=ft[
    ft["table_type"].eq("transition")
    & ft["parent_level"].eq("broad")
    & ft["child_level"].eq("class_name")
]
cal=calibrate_signal_detectability(
    q,biological_unit_col="heldout_unit",
    A_col="A",B_col="B",C_col="C",
    fractions=[0,.05,.1,.2]
)

extra={
    "bootstrap":boot,
    "ident":ident,
    "metric":metric,
    "feature_sensitivity":fs,
    "permutation":perm,
    "calibration":cal,
}
lines=build_interpretation(fit,extra)
assert len(lines)>=5
assert any("selected" in x.lower() for x in lines)

report["full"]="PASS"
report["selected_level"]=fit.selected_level
report["adjusted_feature_counts"]=counts
report["interpretation_lines"]=len(lines)

print(json.dumps(report,indent=2))
