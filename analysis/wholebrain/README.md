Whole-brain manuscript reproduction materials
This directory contains manuscript-specific audit/reproduction materials for the
whole-brain Broad→Class analyses reported in the TransferRes manuscript.
These files are separate from the frozen TransferRes v1.0.0 software release.
They do not modify the v1.0.0 estimator, statistical contract, Zenodo archive, or
the existing Whole-brain Glut Broad→Class 14/14 release gate.
Files
`common_taxonomy_permutation_audit_FINAL.csv`
Frozen common-taxonomy audit table used for Extended Data Fig. 3. It includes
the glutamatergic and GABAergic Broad→Class rows and the downstream hierarchy
transitions.
`broad_class_neuronal_family_summary.csv`
Machine-readable manuscript summary for the glutamatergic and GABAergic
Broad→Class comparator, including raw held-out gain, partial-pooling gain,
shrinkage weight, geometry and reported permutation P value.
`verify_gaba_broad_to_class.py`
Lightweight verification script that checks consistency between the frozen
audit table and the manuscript summary and verifies the GABAergic values
reported in the manuscript.
GABAergic Broad→Class geometry
The GABAergic analysis used 26 animals, 52 perturbations and eight Classes:
06 CTX-CGE GABA
07 CTX-MGE GABA
08 CNU-MGE GABA
09 CNU-LGE GABA
11 CNU-HYa GABA
12 HY GABA
20 MB GABA
28 CB GABA
The frozen manuscript summary contains 6,401 matched
source × perturbation × Class units.
Reported GABAergic Broad→Class results:
raw held-out gain: -5.09%
partial-pooling gain: 0.000%
mean training-derived shrinkage weight: 0.00%
common-taxonomy permutation P: 0.465
The frozen common-taxonomy audit used 100 GABAergic permutations and records
`P_greater = 0.465347`, which rounds to the manuscript value of 0.465.
Run
From the repository root:
```bash
python analysis/wholebrain/verify_gaba_broad_to_class.py
```
A successful run ends with:
```text
OVERALL: PASS
```
Scope
This directory is a manuscript audit/reproduction layer. The verification script
checks the frozen reported outputs against the archived audit table; it does not
reconstruct cell-level perturbational responses or rerun the full
common-taxonomy permutation procedure from raw whole-brain Perturb-seq data.
The core estimator and statistical procedures are implemented in the
`transferres/` package. The whole-brain analysis design and data-processing
details are described in the manuscript Methods.
