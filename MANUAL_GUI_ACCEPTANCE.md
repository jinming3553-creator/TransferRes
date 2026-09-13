# TransferRes Windows GUI Acceptance — manual final check

The automated acceptance script validates the packaged EXE and statistical core.
After it reports `OVERALL: PASS`, do these four visual checks:

1. Double-click `dist\TransferRes.exe`.
2. Confirm the window opens without a console window or error dialog.
3. Click `Load Sample Data` → `Recommended` → `Run Analysis`.
4. Confirm:
   - progress bar advances;
   - Results tab opens;
   - selected resolution is `class_name`;
   - `Open HTML Report` works;
   - `resolution_curve.png` is present.

## Display scaling
Repeat only the window-opening check at:
- 100%
- 125%
- 150%

Pass criterion: all controls remain accessible and text is not unusably clipped.

If all automated tests + these visual checks pass, v0.6.1 RC can be promoted to
the v1.0 packaging stage without changing the statistical estimator.
