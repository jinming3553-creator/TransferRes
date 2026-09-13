# TransferRes

TransferRes is a point-and-click desktop tool for estimating the resolution at
which perturbational responses transfer across biological units.

## Recommended first run

1. Launch TransferRes.
2. Click **Load Sample Data**.
3. Choose **Recommended**.
4. Click **Validate Inputs**.
5. Click **Run Analysis**.
6. Review the Results tab and HTML report.

## Core interpretation

TransferRes distinguishes:
- the discrete resolution that minimizes raw held-out SSE;
- a shrunken finer-level transferable component;
- structural identifiability;
- detectability under the observed sampling geometry.

These quantities answer different questions and are intentionally not collapsed
into a single score.

## Windows release build

Run:

`Build_Release_Windows.bat`

A portable ZIP is created only after the automated Windows acceptance suite
passes.
