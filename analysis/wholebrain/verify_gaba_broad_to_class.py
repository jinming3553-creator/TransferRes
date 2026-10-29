#!/usr/bin/env python3
"""
Verify the frozen whole-brain Broad→Class neuronal-family manuscript outputs.

This is a manuscript audit/reproduction wrapper. It verifies consistency between:
  1) common_taxonomy_permutation_audit_FINAL.csv
  2) broad_class_neuronal_family_summary.csv

It does not reconstruct cell-level responses or rerun the original permutations
from raw whole-brain Perturb-seq data.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
AUDIT = HERE / "common_taxonomy_permutation_audit_FINAL.csv"
SUMMARY = HERE / "broad_class_neuronal_family_summary.csv"


def fail(message: str) -> None:
    print(f"[FAIL] {message}")
    raise SystemExit(1)


def pass_(message: str) -> None:
    print(f"[PASS] {message}")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        fail(f"Missing required file: {path.name}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def i(row: dict[str, str], key: str) -> int:
    return int(row[key])


def assert_close(label: str, observed: float, expected: float, atol: float) -> None:
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=atol):
        fail(f"{label}: observed={observed!r}, expected={expected!r}, atol={atol}")
    pass_(f"{label}: {observed}")


def assert_equal(label: str, observed, expected) -> None:
    if observed != expected:
        fail(f"{label}: observed={observed!r}, expected={expected!r}")
    pass_(f"{label}: {observed}")


def main() -> int:
    audit_rows = read_rows(AUDIT)
    summary_rows = read_rows(SUMMARY)

    audit = {r["transition"]: r for r in audit_rows}
    summary = {r["population"]: r for r in summary_rows}

    required_audit = {
        "Broad→Class (Glutamatergic)",
        "Broad→Class (GABAergic)",
    }
    required_summary = {"Glutamatergic", "GABAergic"}

    missing_audit = required_audit - set(audit)
    missing_summary = required_summary - set(summary)
    if missing_audit:
        fail(f"Missing audit rows: {sorted(missing_audit)}")
    if missing_summary:
        fail(f"Missing summary rows: {sorted(missing_summary)}")

    g = summary["GABAergic"]
    ga = audit["Broad→Class (GABAergic)"]

    # Frozen manuscript geometry.
    assert_equal("GABA animals", i(g, "n_animals"), 26)
    assert_equal("GABA perturbations", i(g, "n_perturbations"), 52)
    assert_equal("GABA Classes", i(g, "n_classes"), 8)
    assert_equal("GABA matched units", i(g, "n_matched_units"), 6401)

    # Frozen manuscript values.
    assert_close("GABA raw held-out gain (%)",
                 f(g, "raw_heldout_gain_percent"), -5.09, 1e-12)
    assert_close("GABA partial-pooling gain (%)",
                 f(g, "partial_pooling_gain_percent"), 0.0, 1e-12)
    assert_close("GABA mean shrinkage weight (%)",
                 f(g, "mean_shrinkage_weight_percent"), 0.0, 1e-12)
    assert_close("GABA reported common-taxonomy P",
                 f(g, "common_taxonomy_p_reported"), 0.465, 5e-4)

    # Consistency with frozen audit table.
    assert_close(
        "GABA summary vs audit partial-pooling gain",
        f(g, "partial_pooling_gain_percent") / 100.0,
        f(ga, "observed_nested_gain"),
        1e-12,
    )
    assert_close(
        "GABA summary vs audit shrinkage weight",
        f(g, "mean_shrinkage_weight_percent") / 100.0,
        f(ga, "mean_nested_weight"),
        1e-12,
    )
    assert_equal("GABA common-taxonomy permutations",
                 i(ga, "n_permutations"), 100)
    assert_close("GABA audit P_greater",
                 f(ga, "P_greater"), 0.465347, 5e-7)

    # Check finite-sample correction compatibility:
    # P = (k + 1) / (B + 1).
    B = i(ga, "n_permutations")
    p = f(ga, "P_greater")
    k = round(p * (B + 1) - 1)
    p_reconstructed = (k + 1) / (B + 1)
    assert_close("GABA finite-sample P reconstruction",
                 p, p_reconstructed, 5e-7)
    pass_(f"GABA implied null realizations >= observed: {k} of {B}")

    # Companion Glut row, included to verify the neuronal-family comparator.
    gl = summary["Glutamatergic"]
    gla = audit["Broad→Class (Glutamatergic)"]
    assert_equal("Glut animals", i(gl, "n_animals"), 26)
    assert_equal("Glut perturbations", i(gl, "n_perturbations"), 52)
    assert_equal("Glut Classes", i(gl, "n_classes"), 7)
    assert_equal("Glut matched units", i(gl, "n_matched_units"), 5915)
    assert_close("Glut raw held-out gain (%)",
                 f(gl, "raw_heldout_gain_percent"), -5.37, 1e-12)

    # Manuscript summary values are rounded; compare using reporting precision.
    assert_close(
        "Glut summary vs audit partial-pooling gain",
        f(gl, "partial_pooling_gain_percent") / 100.0,
        f(gla, "observed_nested_gain"),
        3e-6,
    )
    assert_close(
        "Glut summary vs audit shrinkage weight",
        f(gl, "mean_shrinkage_weight_percent") / 100.0,
        f(gla, "mean_nested_weight"),
        5e-5,
    )
    assert_equal("Glut common-taxonomy permutations",
                 i(gla, "n_permutations"), 200)
    assert_close("Glut audit P_greater",
                 f(gla, "P_greater"), 0.004975, 5e-7)

    print()
    print("OVERALL: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
