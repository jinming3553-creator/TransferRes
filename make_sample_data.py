
from pathlib import Path
import numpy as np
import pandas as pd

out=Path(__file__).resolve().parent/"sample_data"
out.mkdir(exist_ok=True)

rng=np.random.default_rng(4)
rows=[]
X=[]

for unit in [f"mouse{i}" for i in range(1,9)]:
    for perturbation in ["geneA","geneB","geneC"]:
        for cls in ["classA","classB"]:
            rows.append({
                "source":unit,
                "gene_target":perturbation,
                "broad":"all",
                "class_name":cls,
                "n_match":40,
            })
            shift=np.zeros(20)
            shift[0]=1.0 if cls=="classA" else -1.0
            X.append(shift+rng.normal(0,0.25,20))

pd.DataFrame(rows).to_csv(out/"sample_metadata.csv",index=False)
np.save(out/"sample_response.npy",np.asarray(X))

print(out)
