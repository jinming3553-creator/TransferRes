# TransferRes v0.6.0 Release Candidate

## First-time test

1. Unzip.
2. Double-click `Launch_TransferRes.bat`.
3. Click **Load Sample Data**.
4. Choose `Recommended` or `Full`.
5. Click **Validate Inputs**.
6. Click **Run Analysis**.

All three profiles have been exercised by the bundled release QA.

## Generate your own input template

Click **Generate Input Template** in the GUI.

It creates:
- `TransferRes_metadata_template.csv`
- `TransferRes_response_template.csv`
- `TransferRes_INPUT_README.txt`

## Software-assisted interpretation

The Results tab and HTML report now contain conservative interpretation text.

The interpreter is deliberately constrained:
- zero transferable gain is not translated into “no biology”;
- a significant partial-pooling permutation result is not described as a large raw fine-model benefit;
- perturbation bootstrap is not called a biological population CI;
- metric plateaus are not rewritten as exact cross-metric agreement.

## Feature-count safety

If requested sensitivity counts exceed the available feature dimension,
TransferRes automatically caps them to the available dimension and logs the adjustment.

Example:
`500,1000,2000` with 20 available features becomes `20`.

## Build standalone EXE

On Windows, double-click:

`Build_TransferRes_EXE.bat`

Output:
`dist\TransferRes.exe`

## Release QA

Run:

`uv run --with numpy --with pandas --with numba python release_QA.py`

Expected:
- Core: PASS
- Recommended: PASS
- Full: PASS
