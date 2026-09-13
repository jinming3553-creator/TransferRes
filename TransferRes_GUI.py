
from __future__ import annotations

import json
import html
import os
import sys
import traceback
import webbrowser
import threading
import queue
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from interpretation import build_interpretation

from transferres import (
    fit_resolution,
    bootstrap_biological_units,
    diagnose_identifiability,
    diagnose_metric_sensitivity,
    diagnose_feature_count_sensitivity,
    common_taxonomy_permutation_test,
    calibrate_signal_detectability,
)

APP_VERSION = "1.0.0"


def load_response(path: str) -> np.ndarray:
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".npy":
        X = np.load(p)
    elif ext == ".npz":
        z = np.load(p)
        if "response" in z.files:
            X = z["response"]
        elif len(z.files) == 1:
            X = z[z.files[0]]
        else:
            raise ValueError(
                "NPZ contains multiple arrays and no array named 'response'. "
                f"Available arrays: {list(z.files)}"
            )
    elif ext == ".csv":
        X = pd.read_csv(p).to_numpy()
    else:
        raise ValueError("Response file must be .npy, .npz, or .csv")

    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"Response matrix must be 2D; got shape {X.shape}")
    return X


def import_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


class CancelledAnalysis(Exception):
    pass



def run_executable_self_test():
    """
    Headless acceptance test intended to run inside the packaged Windows EXE.

    Returns process exit code:
      0 = PASS
      1 = FAIL
    """
    try:
        sample_root = ROOT / "sample_data"
        meta_path = sample_root / "sample_metadata.csv"
        response_path = sample_root / "sample_response.npy"

        if not meta_path.exists() or not response_path.exists():
            raise RuntimeError("Bundled sample data missing from executable package.")

        meta = pd.read_csv(meta_path)
        X = np.load(response_path)
        meta = meta.rename(columns={
            "source":"biological_unit",
            "gene_target":"perturbation",
        })
        hierarchy=["broad","class_name"]

        # Core
        fit = fit_resolution(
            X, meta,
            hierarchy=hierarchy,
            weight_col="n_match",
            nested_feature_selection=True,
        )
        if fit.selected_level != "class_name":
            raise AssertionError(
                f"Unexpected selected level: {fit.selected_level}"
            )

        # Recommended
        boot = bootstrap_biological_units(
            fit, n_boot=100, seed=1
        )
        ident = diagnose_identifiability(
            meta,
            parent_col="broad",
            child_col="class_name",
            biological_unit_col="biological_unit",
            perturbation_col="perturbation",
        )
        if boot.selection_frequency.empty:
            raise AssertionError("Bootstrap output is empty.")
        if ident.support_frontier.empty:
            raise AssertionError("Identifiability output is empty.")

        # Full diagnostics, intentionally small N for executable QA.
        metric = diagnose_metric_sensitivity(
            X, meta, hierarchy,
            weight_col="n_match",
            primary_level=fit.selected_level,
            n_boot=100,
            seed=1,
        )

        requested=[500,1000,2000]
        p=X.shape[1]
        counts=sorted({min(x,p) for x in requested})
        fs = diagnose_feature_count_sensitivity(
            X, meta, hierarchy,
            feature_counts=counts,
            weight_col="n_match",
            reference_level=fit.selected_level,
        )

        perm = common_taxonomy_permutation_test(
            X, meta,
            parent_col="broad",
            child_col="class_name",
            weight_col="n_match",
            n_permutations=49,
            seed=4,
        )

        ft=fit.fold_table
        q=ft[
            ft["table_type"].eq("transition")
            & ft["parent_level"].eq("broad")
            & ft["child_level"].eq("class_name")
        ].copy()

        cal=calibrate_signal_detectability(
            q,
            biological_unit_col="heldout_unit",
            A_col="A",
            B_col="B",
            C_col="C",
            fractions=[0,.05,.1,.2],
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
        if len(lines) < 5:
            raise AssertionError("Interpretation output unexpectedly short.")

        # Statistical-contract assertions.
        if fit.metadata.get("resolution_selection_rule") != "minimum raw held-out SSE":
            raise AssertionError("Resolution selection rule changed.")
        if fit.metadata.get("partial_pooling_role") != "diagnostic only; not used to choose selected_level":
            raise AssertionError("Partial-pooling role changed.")
        if not boot.metadata.get("bootstrap_unit") == "biological_unit":
            raise AssertionError("Bootstrap unit changed.")
        if not perm.metadata.get("plus_one_correction"):
            raise AssertionError("Permutation +1 correction missing.")
        if "shared across perturbations" not in perm.metadata.get("mapping_scope",""):
            raise AssertionError("Common-taxonomy permutation mapping scope changed.")

        print("TRANSFERRES_SELF_TEST=PASS")
        print(f"VERSION={APP_VERSION}")
        print(f"SELECTED_LEVEL={fit.selected_level}")
        print(f"N_BIOLOGICAL_UNITS={meta.biological_unit.nunique()}")
        print(f"ADJUSTED_FEATURE_COUNTS={counts}")
        print(f"INTERPRETATION_LINES={len(lines)}")
        return 0

    except Exception as e:
        print("TRANSFERRES_SELF_TEST=FAIL")
        print(f"VERSION={APP_VERSION}")
        print(f"ERROR={e}")
        traceback.print_exc()
        return 1



def dataframe_to_html_table(
    df: pd.DataFrame,
    *,
    selected_level: str | None = None,
    level_column: str = "level",
) -> str:
    """
    Render a DataFrame as explicit semantic HTML.

    This avoids relying on browser/text-extractor handling of pandas' default
    table markup and guarantees one <th> per column.

    Long tables are wrapped in a horizontally scrollable container.
    """
    if df is None:
        return ""

    cols=list(df.columns)

    def fmt(v):
        if pd.isna(v):
            return "NA"
        if isinstance(v, (np.floating, float)):
            av=abs(float(v))
            if av != 0 and (av < 1e-4 or av >= 1e5):
                return f"{float(v):.4e}"
            return f"{float(v):.6g}"
        if isinstance(v, (np.integer, int)):
            return str(int(v))
        return str(v)

    out=[]
    out.append('<div class="table-scroll">')
    out.append('<table class="tr-table">')
    out.append("<thead><tr>")
    for c in cols:
        out.append(f"<th>{html.escape(str(c))}</th>")
    out.append("</tr></thead>")
    out.append("<tbody>")

    for _,row in df.iterrows():
        cls=""
        if (
            selected_level is not None
            and level_column in df.columns
            and str(row[level_column]) == str(selected_level)
        ):
            cls=' class="selected-row"'
        out.append(f"<tr{cls}>")
        for c in cols:
            out.append(f"<td>{html.escape(fmt(row[c]))}</td>")
        out.append("</tr>")

    out.append("</tbody></table></div>")
    return "".join(out)


class TransferResApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"TransferRes {APP_VERSION}")
        self.geometry("1160x880")
        self.minsize(1020, 760)

        self.meta_df = None
        self.response = None
        self.last_output = None

        self.worker = None
        self.cancel_event = threading.Event()
        self.ui_queue = queue.Queue()

        self.meta_path = tk.StringVar()
        self.response_path = tk.StringVar()
        self.output_dir = tk.StringVar(
            value=str(Path.home() / "Documents" / "TransferRes_results")
        )

        self.bio_col = tk.StringVar()
        self.pert_col = tk.StringVar()
        self.weight_col = tk.StringVar(value="<none>")
        self.hierarchy_var = tk.StringVar()

        self.feature_count = tk.StringVar(value="")
        self.feature_method = tk.StringVar(value="variance")
        self.nested_var = tk.BooleanVar(value=True)

        self.bootstrap_var = tk.BooleanVar(value=True)
        self.bootstrap_n = tk.StringVar(value="2000")
        self.ident_var = tk.BooleanVar(value=True)

        self.metric_var = tk.BooleanVar(value=False)
        self.metric_boot = tk.StringVar(value="5000")

        self.feature_sens_var = tk.BooleanVar(value=False)
        self.feature_counts_var = tk.StringVar(value="500,1000,2000")

        self.perm_var = tk.BooleanVar(value=False)
        self.perm_n = tk.StringVar(value="1000")
        self.perm_transition = tk.StringVar(value="<last transition>")

        self.cal_var = tk.BooleanVar(value=False)
        self.cal_fractions = tk.StringVar(
            value="0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.3"
        )

        self.run_profile = tk.StringVar(value="Recommended")
        self.progress_var = tk.DoubleVar(value=0.0)
        self.step_var = tk.StringVar(value="Ready")

        self._build_ui()
        self.after(100, self._process_ui_queue)

    # ---------------- UI ----------------
    def _build_ui(self):
        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")
        ttk.Label(
            top,
            text="TransferRes — Transferable Perturbational Resolution",
            font=("Segoe UI", 17, "bold")
        ).pack(anchor="w")
        ttk.Label(
            top,
            text="Leakage-aware biological-unit holdout analysis with inference and diagnostics."
        ).pack(anchor="w", pady=(2, 4))

        self.tabs = ttk.Notebook(self)
        self.tabs.pack(fill="both", expand=True, padx=8, pady=4)

        # Setup tab uses a scrollable canvas so all controls remain accessible
        # on laptops / high-DPI displays / 125-150% Windows scaling.
        self.setup_tab_container = ttk.Frame(self.tabs)
        self.results_tab = ttk.Frame(self.tabs, padding=10)
        self.tabs.add(self.setup_tab_container, text="Setup & Run")
        self.tabs.add(self.results_tab, text="Results")

        self.setup_canvas = tk.Canvas(
            self.setup_tab_container,
            highlightthickness=0,
            borderwidth=0
        )
        self.setup_scrollbar = ttk.Scrollbar(
            self.setup_tab_container,
            orient="vertical",
            command=self.setup_canvas.yview
        )
        self.setup_canvas.configure(yscrollcommand=self.setup_scrollbar.set)

        self.setup_scrollbar.pack(side="right", fill="y")
        self.setup_canvas.pack(side="left", fill="both", expand=True)

        self.setup_tab = ttk.Frame(self.setup_canvas, padding=10)
        self._setup_window = self.setup_canvas.create_window(
            (0, 0),
            window=self.setup_tab,
            anchor="nw"
        )

        def _sync_setup_scrollregion(event=None):
            self.setup_canvas.configure(scrollregion=self.setup_canvas.bbox("all"))

        def _sync_setup_width(event):
            self.setup_canvas.itemconfigure(self._setup_window, width=event.width)

        self.setup_tab.bind("<Configure>", _sync_setup_scrollregion)
        self.setup_canvas.bind("<Configure>", _sync_setup_width)

        # Mouse wheel support while pointer is over the setup page.
        def _on_setup_mousewheel(event):
            if event.delta:
                self.setup_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _bind_setup_wheel(event):
            self.setup_canvas.bind_all("<MouseWheel>", _on_setup_mousewheel)

        def _unbind_setup_wheel(event):
            self.setup_canvas.unbind_all("<MouseWheel>")

        self.setup_canvas.bind("<Enter>", _bind_setup_wheel)
        self.setup_canvas.bind("<Leave>", _unbind_setup_wheel)

        self._build_setup()
        self._build_results()

    def _build_setup(self):
        main = self.setup_tab

        files = ttk.LabelFrame(main, text="1. Input files", padding=10)
        files.pack(fill="x", pady=5)

        btnrow = ttk.Frame(files)
        btnrow.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        ttk.Button(btnrow, text="Load Sample Data", command=self._load_sample_data).pack(side="left")
        ttk.Button(btnrow, text="Generate Input Template", command=self._generate_template).pack(side="left", padx=8)
        ttk.Button(btnrow, text="Clear Inputs", command=self._clear_inputs).pack(side="left", padx=8)

        self._file_row(files, 1, "Metadata CSV", self.meta_path, self._browse_meta)
        self._file_row(files, 2, "Response matrix (.npy/.npz/.csv)", self.response_path, self._browse_response)
        self._file_row(files, 3, "Output folder", self.output_dir, self._browse_output)

        mapping = ttk.LabelFrame(main, text="2. Metadata mapping", padding=10)
        mapping.pack(fill="x", pady=5)

        ttk.Label(mapping, text="Biological unit").grid(row=0, column=0, sticky="w")
        self.bio_combo = ttk.Combobox(mapping, textvariable=self.bio_col, state="readonly", width=28)
        self.bio_combo.grid(row=0, column=1, sticky="w", padx=8)

        ttk.Label(mapping, text="Perturbation").grid(row=0, column=2, sticky="w")
        self.pert_combo = ttk.Combobox(mapping, textvariable=self.pert_col, state="readonly", width=28)
        self.pert_combo.grid(row=0, column=3, sticky="w", padx=8)

        ttk.Label(mapping, text="Weight / depth").grid(row=1, column=0, sticky="w", pady=6)
        self.weight_combo = ttk.Combobox(mapping, textvariable=self.weight_col, state="readonly", width=28)
        self.weight_combo.grid(row=1, column=1, sticky="w", padx=8)

        ttk.Label(mapping, text="Hierarchy coarse → fine").grid(row=2, column=0, sticky="w", pady=6)
        ttk.Entry(mapping, textvariable=self.hierarchy_var, width=85).grid(
            row=2, column=1, columnspan=3, sticky="ew", padx=8
        )
        ttk.Label(
            mapping,
            text="Comma-separated metadata columns, e.g. broad,class,subclass,supertype,cluster",
            foreground="#555"
        ).grid(row=3, column=1, columnspan=3, sticky="w", padx=8)
        mapping.columnconfigure(3, weight=1)

        core = ttk.LabelFrame(main, text="3. Core analysis", padding=10)
        core.pack(fill="x", pady=5)

        ttk.Label(core, text="Run profile").grid(row=0, column=0, sticky="w")
        profile = ttk.Combobox(
            core,
            textvariable=self.run_profile,
            values=["Core", "Recommended", "Full"],
            state="readonly",
            width=18
        )
        profile.grid(row=0, column=1, sticky="w", padx=8)
        profile.bind("<<ComboboxSelected>>", lambda e: self._apply_profile())

        ttk.Label(core, text="Feature count").grid(row=1, column=0, sticky="w", pady=5)
        ttk.Entry(core, textvariable=self.feature_count, width=12).grid(
            row=1, column=1, sticky="w", padx=8
        )
        ttk.Label(core, text="blank = all").grid(row=1, column=2, sticky="w")

        ttk.Label(core, text="Feature ranking").grid(row=2, column=0, sticky="w", pady=5)
        ttk.Combobox(
            core,
            textvariable=self.feature_method,
            values=["variance", "detection"],
            state="readonly",
            width=16
        ).grid(row=2, column=1, sticky="w", padx=8)

        ttk.Checkbutton(
            core,
            text="Nested feature selection (recommended; leakage-free)",
            variable=self.nested_var
        ).grid(row=3, column=0, columnspan=3, sticky="w", pady=3)

        inf = ttk.LabelFrame(main, text="4. Inference & diagnostics", padding=10)
        inf.pack(fill="x", pady=5)

        ttk.Checkbutton(inf, text="Biological-unit bootstrap", variable=self.bootstrap_var).grid(
            row=0, column=0, sticky="w"
        )
        ttk.Label(inf, text="N").grid(row=0, column=1, sticky="e")
        ttk.Entry(inf, textvariable=self.bootstrap_n, width=9).grid(
            row=0, column=2, sticky="w", padx=5
        )

        ttk.Checkbutton(inf, text="Structural identifiability", variable=self.ident_var).grid(
            row=1, column=0, sticky="w", pady=3
        )

        ttk.Checkbutton(inf, text="Metric sensitivity", variable=self.metric_var).grid(
            row=2, column=0, sticky="w", pady=3
        )
        ttk.Label(inf, text="Perturbation bootstrap N").grid(row=2, column=1, sticky="e")
        ttk.Entry(inf, textvariable=self.metric_boot, width=9).grid(
            row=2, column=2, sticky="w", padx=5
        )

        ttk.Checkbutton(inf, text="Feature-count sensitivity", variable=self.feature_sens_var).grid(
            row=3, column=0, sticky="w", pady=3
        )
        ttk.Label(inf, text="Counts").grid(row=3, column=1, sticky="e")
        ttk.Entry(inf, textvariable=self.feature_counts_var, width=22).grid(
            row=3, column=2, sticky="w", padx=5
        )

        ttk.Checkbutton(inf, text="Common-taxonomy permutation", variable=self.perm_var).grid(
            row=4, column=0, sticky="w", pady=3
        )
        ttk.Label(inf, text="N").grid(row=4, column=1, sticky="e")
        ttk.Entry(inf, textvariable=self.perm_n, width=9).grid(
            row=4, column=2, sticky="w", padx=5
        )
        ttk.Label(inf, text="Transition").grid(row=4, column=3, sticky="e", padx=(15, 2))
        self.perm_combo = ttk.Combobox(
            inf, textvariable=self.perm_transition, state="readonly", width=28
        )
        self.perm_combo.grid(row=4, column=4, sticky="w")

        ttk.Checkbutton(inf, text="Signal-injection calibration", variable=self.cal_var).grid(
            row=5, column=0, sticky="w", pady=3
        )
        ttk.Label(inf, text="Fractions").grid(row=5, column=1, sticky="e")
        ttk.Entry(inf, textvariable=self.cal_fractions, width=45).grid(
            row=5, column=2, columnspan=3, sticky="w", padx=5
        )

        controls = ttk.Frame(main)
        controls.pack(fill="x", pady=(8, 4))
        self.run_btn = ttk.Button(controls, text="Run Analysis", command=self._start_run)
        self.run_btn.pack(side="left")
        self.cancel_btn = ttk.Button(controls, text="Cancel", command=self._cancel_run, state="disabled")
        self.cancel_btn.pack(side="left", padx=8)
        ttk.Button(controls, text="Validate Inputs", command=self._validate_only).pack(side="left", padx=8)
        ttk.Button(controls, text="Open Output Folder", command=self._open_output).pack(side="left")
        ttk.Button(controls, text="About", command=self._show_about).pack(side="left", padx=8)

        progress = ttk.Frame(main)
        progress.pack(fill="x", pady=(2, 4))
        ttk.Progressbar(progress, variable=self.progress_var, maximum=100).pack(
            side="left", fill="x", expand=True
        )
        ttk.Label(progress, textvariable=self.step_var, width=36).pack(side="left", padx=8)

        sf = ttk.LabelFrame(main, text="Status / warnings", padding=8)
        sf.pack(fill="both", expand=True, pady=5)
        self.status = tk.Text(sf, height=10, wrap="word")
        self.status.pack(fill="both", expand=True)
        self._log("Ready. One metadata row must represent one aggregated analysis unit, not one raw cell.\n")
        self.after_idle(
            lambda: self.setup_canvas.configure(scrollregion=self.setup_canvas.bbox("all"))
        )

    def _build_results(self):
        main = self.results_tab

        summary = ttk.LabelFrame(main, text="Analysis summary", padding=10)
        summary.pack(fill="x", pady=5)

        self.result_selected = tk.StringVar(value="No analysis run yet.")
        ttk.Label(summary, text="Selected resolution:", font=("Segoe UI", 11, "bold")).grid(
            row=0, column=0, sticky="w"
        )
        ttk.Label(summary, textvariable=self.result_selected, font=("Segoe UI", 11)).grid(
            row=0, column=1, sticky="w", padx=8
        )

        self.result_path = tk.StringVar(value="")
        ttk.Label(summary, text="Result folder:").grid(row=1, column=0, sticky="w", pady=3)
        ttk.Label(summary, textvariable=self.result_path).grid(row=1, column=1, sticky="w", padx=8)

        actions = ttk.Frame(summary)
        actions.grid(row=2, column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Button(actions, text="Open Results Folder", command=self._open_last_result).pack(side="left")
        ttk.Button(actions, text="Open HTML Report", command=self._open_html_report).pack(
            side="left", padx=8
        )
        ttk.Button(actions, text="Open Error Log", command=self._open_error_log).pack(
            side="left", padx=8
        )

        self.interpretation_box = tk.Text(summary, height=6, wrap="word")
        self.interpretation_box.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(10,0))

        pane = ttk.Panedwindow(main, orient="horizontal")
        pane.pack(fill="both", expand=True, pady=5)

        left = ttk.LabelFrame(pane, text="Key tables", padding=6)
        right = ttk.LabelFrame(pane, text="Plots", padding=6)
        pane.add(left, weight=3)
        pane.add(right, weight=2)

        self.result_text = tk.Text(left, wrap="none")
        self.result_text.pack(fill="both", expand=True)

        self.plot_list = tk.Listbox(right)
        self.plot_list.pack(fill="both", expand=True)
        ttk.Button(right, text="Open selected plot", command=self._open_selected_plot).pack(
            fill="x", pady=(6, 0)
        )

    # ---------------- profiles ----------------
    def _apply_profile(self):
        p = self.run_profile.get()
        if p == "Core":
            self.bootstrap_var.set(False)
            self.ident_var.set(False)
            self.metric_var.set(False)
            self.feature_sens_var.set(False)
            self.perm_var.set(False)
            self.cal_var.set(False)
        elif p == "Recommended":
            self.bootstrap_var.set(True)
            self.ident_var.set(True)
            self.metric_var.set(False)
            self.feature_sens_var.set(False)
            self.perm_var.set(False)
            self.cal_var.set(False)
        else:
            self.bootstrap_var.set(True)
            self.ident_var.set(True)
            self.metric_var.set(True)
            self.feature_sens_var.set(True)
            self.perm_var.set(True)
            self.cal_var.set(True)

    # ---------------- sample ----------------
    def _load_sample_data(self):
        meta = ROOT / "sample_data" / "sample_metadata.csv"
        resp = ROOT / "sample_data" / "sample_response.npy"
        if not meta.exists() or not resp.exists():
            messagebox.showerror("Sample data", "Bundled sample data not found.")
            return

        self.meta_path.set(str(meta))
        self.response_path.set(str(resp))
        self.meta_df = pd.read_csv(meta)
        self.response = np.load(resp)

        cols = list(self.meta_df.columns)
        self.bio_combo["values"] = cols
        self.pert_combo["values"] = cols
        self.weight_combo["values"] = ["<none>"] + cols

        self.bio_col.set("source")
        self.pert_col.set("gene_target")
        self.weight_col.set("n_match")
        self.hierarchy_var.set("broad,class_name")
        self.feature_count.set("")
        self.run_profile.set("Recommended")
        self._apply_profile()

        self._log("Loaded bundled sample data.\n")


    def _generate_template(self):
        dest=filedialog.askdirectory(title="Choose folder for TransferRes input template")
        if not dest:
            return
        d=Path(dest)
        meta=pd.DataFrame([
            {
                "biological_unit":"unit1",
                "perturbation":"perturbationA",
                "broad":"all",
                "class_name":"classA",
                "n_match":20,
            },
            {
                "biological_unit":"unit2",
                "perturbation":"perturbationA",
                "broad":"all",
                "class_name":"classA",
                "n_match":20,
            },
        ])
        meta.to_csv(d/"TransferRes_metadata_template.csv",index=False)

        response_example=pd.DataFrame({
            "feature_1":[0.1,0.2],
            "feature_2":[-0.1,0.0],
            "feature_3":[0.3,0.4],
        })
        response_example.to_csv(d/"TransferRes_response_template.csv",index=False)

        (d/"TransferRes_INPUT_README.txt").write_text(
            "TransferRes input template\n\n"
            "1. Metadata and response files must have exactly the same number of rows.\n"
            "2. Each row must represent one aggregated analysis unit, not one raw cell.\n"
            "3. biological_unit must identify independent biological replicates.\n"
            "4. perturbation identifies the perturbation condition.\n"
            "5. hierarchy columns must be nested from coarse to fine.\n"
            "6. response columns are numeric response features.\n"
            "7. Do not duplicate biological_unit × perturbation × hierarchy rows.\n",
            encoding="utf-8",
        )
        messagebox.showinfo("Template created",f"Input templates saved to:\n{d}")

    def _clear_inputs(self):
        self.meta_path.set("")
        self.response_path.set("")
        self.meta_df = None
        self.response = None
        self.bio_col.set("")
        self.pert_col.set("")
        self.weight_col.set("<none>")
        self.hierarchy_var.set("")
        self._log("Inputs cleared.\n")

    # ---------------- input helpers ----------------
    def _file_row(self, parent, row, label, var, command):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=3)
        ttk.Entry(parent, textvariable=var, width=85).grid(
            row=row, column=1, sticky="ew", padx=8
        )
        ttk.Button(parent, text="Browse…", command=command).grid(row=row, column=2)
        parent.columnconfigure(1, weight=1)

    def _browse_meta(self):
        p = filedialog.askopenfilename(
            title="Choose metadata CSV",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")]
        )
        if not p:
            return
        self.meta_path.set(p)
        try:
            self.meta_df = pd.read_csv(p)
            cols = list(self.meta_df.columns)
            self.bio_combo["values"] = cols
            self.pert_combo["values"] = cols
            self.weight_combo["values"] = ["<none>"] + cols
            for cand in ["biological_unit", "source", "donor", "replicate", "mouse"]:
                if cand in cols:
                    self.bio_col.set(cand)
                    break
            for cand in ["perturbation", "gene_target", "drug", "guide"]:
                if cand in cols:
                    self.pert_col.set(cand)
                    break
            for cand in ["n_match", "weight", "depth", "n_cells"]:
                if cand in cols:
                    self.weight_col.set(cand)
                    break
            self._log(
                f"Loaded metadata: {self.meta_df.shape[0]} rows × "
                f"{self.meta_df.shape[1]} columns\n"
            )
        except Exception as e:
            messagebox.showerror("Metadata error", str(e))

    def _browse_response(self):
        p = filedialog.askopenfilename(
            title="Choose response matrix",
            filetypes=[("NumPy", "*.npy *.npz"), ("CSV", "*.csv"), ("All files", "*.*")]
        )
        if not p:
            return
        self.response_path.set(p)
        try:
            self.response = load_response(p)
            self._log(f"Loaded response matrix: shape={self.response.shape}\n")
        except Exception as e:
            messagebox.showerror("Response error", str(e))

    def _browse_output(self):
        p = filedialog.askdirectory(title="Choose output folder")
        if p:
            self.output_dir.set(p)

    def _parse_config(self):
        if self.meta_df is None:
            if not self.meta_path.get():
                raise ValueError("Choose metadata CSV.")
            self.meta_df = pd.read_csv(self.meta_path.get())

        if self.response is None:
            if not self.response_path.get():
                raise ValueError("Choose response matrix.")
            self.response = load_response(self.response_path.get())

        if len(self.meta_df) != len(self.response):
            raise ValueError(
                f"Metadata rows ({len(self.meta_df)}) != response rows ({len(self.response)})."
            )

        bio = self.bio_col.get().strip()
        pert = self.pert_col.get().strip()
        if not bio or not pert:
            raise ValueError("Select biological-unit and perturbation columns.")

        hierarchy = [
            x.strip() for x in self.hierarchy_var.get().split(",") if x.strip()
        ]
        if not hierarchy:
            raise ValueError("Enter hierarchy columns.")

        missing = [c for c in [bio, pert] + hierarchy if c not in self.meta_df.columns]
        if missing:
            raise ValueError(f"Columns not found: {missing}")

        w = self.weight_col.get().strip()
        weight = None if w in ("", "<none>") else w

        meta = self.meta_df.copy()
        rename = {}
        if bio != "biological_unit":
            rename[bio] = "biological_unit"
        if pert != "perturbation":
            rename[pert] = "perturbation"
        meta = meta.rename(columns=rename)
        hierarchy = [rename.get(c, c) for c in hierarchy]
        weight = rename.get(weight, weight) if weight else None

        if meta["biological_unit"].nunique() < 2:
            raise ValueError("At least two independent biological units are required.")

        dup_cols = ["biological_unit", "perturbation"] + hierarchy
        if meta.duplicated(subset=dup_cols).any():
            raise ValueError(
                "Duplicate biological_unit × perturbation × hierarchy rows detected. "
                "Aggregate raw cells before analysis."
            )

        fc = self.feature_count.get().strip()
        n_features = None if fc == "" else int(fc)
        if n_features is not None and n_features <= 0:
            raise ValueError("Feature count must be positive.")

        trans = [f"{a} → {b}" for a, b in zip(hierarchy[:-1], hierarchy[1:])]
        self.perm_combo["values"] = ["<last transition>"] + trans
        if self.perm_transition.get() not in self.perm_combo["values"]:
            self.perm_transition.set("<last transition>")

        return meta, hierarchy, weight, n_features

    def _validate_only(self):
        try:
            meta, hierarchy, weight, n_features = self._parse_config()
            self._log("\nVALIDATION PASSED\n")
            self._log(f"Biological units: {meta.biological_unit.nunique()}\n")
            self._log(f"Perturbations: {meta.perturbation.nunique()}\n")
            self._log(f"Hierarchy: {' → '.join(hierarchy)}\n")
            if meta.biological_unit.nunique() < 5:
                self._log(
                    "WARNING: <5 biological units; bootstrap intervals are stability diagnostics.\n"
                )
            messagebox.showinfo("Validation", "Inputs passed validation.")
        except Exception as e:
            self._log("\nVALIDATION FAILED\n" + str(e) + "\n")
            messagebox.showerror("Validation failed", str(e))

    # ---------------- threading ----------------
    def _start_run(self):
        try:
            config = self._parse_config()
        except Exception as e:
            messagebox.showerror("Validation failed", str(e))
            return

        if self.worker and self.worker.is_alive():
            return

        self.cancel_event.clear()
        self.run_btn.config(state="disabled")
        self.cancel_btn.config(state="normal")
        self.progress_var.set(0)
        self.step_var.set("Starting…")

        self.worker = threading.Thread(
            target=self._run_worker,
            args=(config,),
            daemon=True
        )
        self.worker.start()

    def _cancel_run(self):
        self.cancel_event.set()
        self.step_var.set("Cancel requested…")
        self._log(
            "Cancel requested. TransferRes will stop after the current calculation step finishes.\n"
        )

    def _check_cancel(self):
        if self.cancel_event.is_set():
            raise CancelledAnalysis("Analysis cancelled by user.")

    def _emit(self, kind, payload=None):
        self.ui_queue.put((kind, payload))

    def _process_ui_queue(self):
        try:
            while True:
                kind, payload = self.ui_queue.get_nowait()
                if kind == "log":
                    self._log(payload)
                elif kind == "progress":
                    val, label = payload
                    self.progress_var.set(val)
                    self.step_var.set(label)
                elif kind == "done":
                    fit, outdir, extra, plots = payload
                    self.last_output = outdir
                    self._populate_results(fit, outdir, extra, plots)
                    self.run_btn.config(state="normal")
                    self.cancel_btn.config(state="disabled")
                    self.progress_var.set(100)
                    self.step_var.set("Done")
                    self.tabs.select(self.results_tab)
                    messagebox.showinfo(
                        "Analysis complete",
                        f"Selected level: {fit.selected_level}\n\nResults:\n{outdir}"
                    )
                elif kind == "cancelled":
                    self.run_btn.config(state="normal")
                    self.cancel_btn.config(state="disabled")
                    self.step_var.set("Cancelled")
                    self._log("Analysis cancelled.\n")
                elif kind == "error":
                    msg, outdir = payload
                    self.run_btn.config(state="normal")
                    self.cancel_btn.config(state="disabled")
                    self.step_var.set("Error")
                    self.last_output = outdir
                    messagebox.showerror("Analysis error", msg)
        except queue.Empty:
            pass
        self.after(100, self._process_ui_queue)

    # ---------------- run worker ----------------
    def _run_worker(self, config):
        meta, hierarchy, weight, n_features = config

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = Path(self.output_dir.get()) / f"TransferRes_{stamp}"
        outdir.mkdir(parents=True, exist_ok=True)

        try:
            self._emit("log", "\n=== RUN STARTED ===\n")
            self._emit("progress", (5, "Core estimator"))

            fit = fit_resolution(
                self.response,
                meta,
                hierarchy=hierarchy,
                weight_col=weight,
                n_features=n_features,
                feature_method=self.feature_method.get(),
                nested_feature_selection=self.nested_var.get(),
            )
            self._check_cancel()

            fit.level_table.to_csv(outdir / "level_table.csv", index=False)
            fit.transition_table.to_csv(outdir / "transition_table.csv", index=False)
            fit.fold_table.to_csv(outdir / "fold_table.csv", index=False)

            self._emit("log", f"Selected transferable resolution: {fit.selected_level}\n")
            extra = {}

            steps = []
            if self.bootstrap_var.get():
                steps.append("bootstrap")
            if self.ident_var.get() and len(hierarchy) >= 2:
                steps.append("ident")
            if self.metric_var.get():
                steps.append("metric")
            if self.feature_sens_var.get():
                steps.append("feature")
            if self.perm_var.get():
                steps.append("perm")
            if self.cal_var.get():
                steps.append("cal")

            total = max(1, len(steps) + 2)
            completed = 1

            def prog(label):
                nonlocal completed
                completed += 1
                pct = min(90, 5 + 85 * completed / total)
                self._emit("progress", (pct, label))

            if "bootstrap" in steps:
                self._check_cancel()
                n = int(self.bootstrap_n.get())
                prog(f"Bootstrap ({n})")
                b = bootstrap_biological_units(fit, n_boot=n, seed=2026)
                b.raw_level_table.to_csv(outdir / "bootstrap_raw_levels.csv", index=False)
                b.transition_table.to_csv(outdir / "bootstrap_transitions.csv", index=False)
                b.selection_frequency.to_csv(outdir / "bootstrap_selection_frequency.csv", index=False)
                extra["bootstrap"] = b
                if b.metadata["low_unit_count_warning"]:
                    self._emit(
                        "log",
                        "WARNING: low unit count; bootstrap is a stability diagnostic.\n"
                    )

            if "ident" in steps:
                self._check_cancel()
                parent, child = hierarchy[-2], hierarchy[-1]
                prog(f"Identifiability {parent} → {child}")
                i = diagnose_identifiability(
                    meta,
                    parent_col=parent,
                    child_col=child,
                    biological_unit_col="biological_unit",
                    perturbation_col="perturbation",
                )
                i.summary_table.to_csv(outdir / "identifiability_summary.csv", index=False)
                i.support_frontier.to_csv(outdir / "identifiability_support_frontier.csv", index=False)
                i.perturbation_table.to_csv(outdir / "identifiability_by_perturbation.csv", index=False)
                i.stratum_table.to_csv(outdir / "identifiability_strata.csv", index=False)
                extra["ident"] = i

            if "metric" in steps:
                self._check_cancel()
                n = int(self.metric_boot.get())
                prog(f"Metric sensitivity ({n})")
                m = diagnose_metric_sensitivity(
                    self.response,
                    meta,
                    hierarchy,
                    weight_col=weight,
                    n_features=n_features,
                    feature_method=self.feature_method.get(),
                    nested_feature_selection=self.nested_var.get(),
                    primary_level=fit.selected_level,
                    n_boot=n,
                    seed=20260913,
                )
                m.bootstrap_table.to_csv(outdir / "metric_sensitivity.csv", index=False)
                m.primary_contrast_table.to_csv(outdir / "metric_contrasts.csv", index=False)
                m.interpretation_table.to_csv(outdir / "metric_interpretation.csv", index=False)
                extra["metric"] = m

            if "feature" in steps:
                self._check_cancel()
                requested_counts = sorted({
                    int(x.strip())
                    for x in self.feature_counts_var.get().split(",")
                    if x.strip()
                })
                p_available=int(self.response.shape[1])
                counts=sorted({min(x,p_available) for x in requested_counts if x>0})
                if counts != requested_counts:
                    self._emit(
                        "log",
                        f"Feature-count sensitivity adjusted to available dimension {p_available}: "
                        f"{requested_counts} → {counts}\n"
                    )
                prog(f"Feature-count sensitivity {counts}")
                fs = diagnose_feature_count_sensitivity(
                    self.response,
                    meta,
                    hierarchy,
                    feature_counts=counts,
                    weight_col=weight,
                    feature_method=self.feature_method.get(),
                    nested_feature_selection=self.nested_var.get(),
                    reference_level=fit.selected_level,
                )
                fs.level_table.to_csv(outdir / "feature_count_levels.csv", index=False)
                fs.selection_table.to_csv(outdir / "feature_count_selection.csv", index=False)
                fs.interpretation_table.to_csv(outdir / "feature_count_interpretation.csv", index=False)
                extra["feature_sensitivity"] = fs

            if "perm" in steps:
                self._check_cancel()
                if len(hierarchy) < 2:
                    raise ValueError("Permutation requires at least two hierarchy levels.")
                choice = self.perm_transition.get()
                if choice == "<last transition>" or not choice:
                    parent, child = hierarchy[-2], hierarchy[-1]
                else:
                    parent, child = [x.strip() for x in choice.split("→")]
                n = int(self.perm_n.get())
                prog(f"Permutation {parent} → {child} ({n})")
                pm = common_taxonomy_permutation_test(
                    self.response,
                    meta,
                    parent_col=parent,
                    child_col=child,
                    weight_col=weight,
                    n_features=n_features,
                    feature_method=self.feature_method.get(),
                    nested_feature_selection=self.nested_var.get(),
                    n_permutations=n,
                    seed=20260913,
                )
                pm.null_table.to_csv(outdir / "permutation_null.csv", index=False)
                pd.DataFrame([{
                    "parent": parent,
                    "child": child,
                    "observed_raw_gain": pm.observed_raw_incremental_gain,
                    "observed_partial_gain": pm.observed_partial_pooling_gain,
                    "p_raw_greater": pm.p_raw_greater,
                    "p_raw_two_sided": pm.p_raw_two_sided,
                    "p_partial_greater": pm.p_partial_greater,
                    "p_partial_two_sided": pm.p_partial_two_sided,
                }]).to_csv(outdir / "permutation_summary.csv", index=False)
                extra["permutation"] = pm

            if "cal" in steps:
                self._check_cancel()
                if len(fit.transition_table) == 0:
                    raise ValueError("Signal calibration requires at least one hierarchy transition.")
                parent, child = hierarchy[-2], hierarchy[-1]
                prog(f"Signal calibration {parent} → {child}")
                ft = fit.fold_table
                q = ft[
                    ft["table_type"].eq("transition")
                    & ft["parent_level"].eq(parent)
                    & ft["child_level"].eq(child)
                ].copy()
                fractions = [
                    float(x.strip())
                    for x in self.cal_fractions.get().split(",")
                    if x.strip()
                ]
                cal = calibrate_signal_detectability(
                    q,
                    biological_unit_col="heldout_unit",
                    A_col="A",
                    B_col="B",
                    C_col="C",
                    fractions=fractions,
                )
                cal.calibration_table.to_csv(outdir / "signal_calibration.csv", index=False)
                extra["calibration"] = cal

            self._check_cancel()
            self._emit("progress", (92, "Exporting report"))

            manifest = {
                "software": {"name": "TransferRes", "version": APP_VERSION},
                "timestamp": datetime.now().isoformat(),
                "input": {
                    "metadata_file": str(Path(self.meta_path.get()).resolve()),
                    "response_file": str(Path(self.response_path.get()).resolve()),
                    "n_rows": int(len(meta)),
                    "n_features_available": int(self.response.shape[1]),
                    "n_biological_units": int(meta.biological_unit.nunique()),
                    "n_perturbations": int(meta.perturbation.nunique()),
                },
                "analysis": {
                    "hierarchy": hierarchy,
                    "weight_col": weight,
                    "feature_count": n_features,
                    "feature_method": self.feature_method.get(),
                    "nested_feature_selection": bool(self.nested_var.get()),
                    "selected_level": fit.selected_level,
                    "selection_rule": "minimum raw held-out SSE",
                    "run_profile": self.run_profile.get(),
                },
                "diagnostics": {
                    "bootstrap": bool(self.bootstrap_var.get()),
                    "identifiability": bool(self.ident_var.get()),
                    "metric_sensitivity": bool(self.metric_var.get()),
                    "feature_count_sensitivity": bool(self.feature_sens_var.get()),
                    "common_taxonomy_permutation": bool(self.perm_var.get()),
                    "signal_calibration": bool(self.cal_var.get()),
                },
                "warnings": {
                    "transductive_preprocessing": not bool(self.nested_var.get()),
                    "low_biological_unit_count": bool(meta.biological_unit.nunique() < 5),
                    "raw_cells_not_supported": True,
                },
                "interpretation_contract": {
                    "zero_gain_does_not_mean_no_biology": True,
                    "partial_pooling_is_diagnostic_not_resolution_selection": True,
                    "perturbation_bootstrap_is_not_biological_population_CI": True,
                },
            }

            (outdir / "analysis_manifest.json").write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            interpretation_lines = build_interpretation(fit, extra)
            manifest["software_interpretation"] = interpretation_lines
            (outdir/"interpretation.txt").write_text("\n".join(f"- {x}" for x in interpretation_lines)+"\n", encoding="utf-8")

            plots = self._make_plots(outdir, fit, extra)
            self._write_html(outdir, fit, extra, plots, manifest)

            self._emit("progress", (100, "Done"))
            self._emit("done", (fit, outdir, extra, plots))

        except CancelledAnalysis:
            self._write_cancel_log(outdir)
            self._emit("cancelled")
        except Exception as e:
            tb = traceback.format_exc()
            self._write_error_log(outdir, e, tb)
            self._emit("log", "\nERROR\n" + str(e) + "\n")
            self._emit("error", (str(e), outdir))

    # ---------------- logging ----------------
    def _write_error_log(self, outdir, exc, tb):
        p = outdir / "error_log.txt"
        p.write_text(
            f"TransferRes {APP_VERSION}\n"
            f"Time: {datetime.now().isoformat()}\n\n"
            f"ERROR:\n{exc}\n\nTRACEBACK:\n{tb}\n",
            encoding="utf-8",
        )

    def _write_cancel_log(self, outdir):
        (outdir / "cancelled.txt").write_text(
            f"TransferRes {APP_VERSION}\n"
            f"Analysis cancelled by user at {datetime.now().isoformat()}.\n",
            encoding="utf-8",
        )

    # ---------------- plots/report ----------------
    def _make_plots(self, outdir, fit, extra):
        plt = import_matplotlib()
        plots = []

        z = fit.level_table
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(np.arange(len(z)), 100 * z["gain_vs_coarsest"].to_numpy(), marker="o")
        ax.axhline(0, linewidth=1)
        ax.set_xticks(np.arange(len(z)), z["level"].tolist(), rotation=25, ha="right")
        ax.set_ylabel("Held-out gain vs coarsest (%)")
        ax.set_title("Transferable resolution curve")
        fig.tight_layout()
        p = outdir / "resolution_curve.png"
        fig.savefig(p, dpi=180, bbox_inches="tight")
        plt.close(fig)
        plots.append(p)

        if len(fit.transition_table):
            t = fit.transition_table
            x = np.arange(len(t))
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(x, 100 * t["raw_incremental_gain"], marker="o", label="Raw fine model")
            ax.plot(x, 100 * t["partial_pooling_gain"], marker="o", label="Partial pooling")
            ax.axhline(0, linewidth=1)
            labels = [f"{a}→{b}" for a, b in zip(t.parent_level, t.child_level)]
            ax.set_xticks(x, labels, rotation=25, ha="right")
            ax.set_ylabel("Incremental gain (%)")
            ax.set_title("Raw vs shrunken transferable component")
            ax.legend(frameon=False)
            fig.tight_layout()
            p = outdir / "transition_gains.png"
            fig.savefig(p, dpi=180, bbox_inches="tight")
            plt.close(fig)
            plots.append(p)

        if extra.get("metric") is not None:
            m = extra["metric"].bootstrap_table
            fig, ax = plt.subplots(figsize=(8, 4.8))
            for metric in ["Euclidean", "Cosine", "Pearson"]:
                q = m[m.metric.eq(metric)]
                ax.plot(np.arange(len(q)), q.mean_score.to_numpy(), marker="o", label=metric)
            q0 = m[m.metric.eq("Euclidean")]
            ax.set_xticks(np.arange(len(q0)), q0.level.tolist(), rotation=25, ha="right")
            ax.set_ylabel("Held-out score")
            ax.set_title("Metric sensitivity")
            ax.legend(frameon=False)
            fig.tight_layout()
            p = outdir / "metric_sensitivity.png"
            fig.savefig(p, dpi=180, bbox_inches="tight")
            plt.close(fig)
            plots.append(p)

        if extra.get("feature_sensitivity") is not None:
            fs = extra["feature_sensitivity"].level_table
            fig, ax = plt.subplots(figsize=(8, 4.8))
            lastg = None
            for nf, g in fs.groupby("n_features"):
                lastg = g
                ax.plot(np.arange(len(g)), 100 * g.gain_vs_coarsest.to_numpy(), marker="o", label=f"{nf} features")
            if lastg is not None:
                ax.set_xticks(np.arange(len(lastg)), lastg.level.tolist(), rotation=25, ha="right")
            ax.axhline(0, linewidth=1)
            ax.set_ylabel("Held-out gain vs coarsest (%)")
            ax.set_title("Feature-count sensitivity")
            ax.legend(frameon=False)
            fig.tight_layout()
            p = outdir / "feature_count_sensitivity.png"
            fig.savefig(p, dpi=180, bbox_inches="tight")
            plt.close(fig)
            plots.append(p)

        if extra.get("calibration") is not None:
            c = extra["calibration"].calibration_table
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(c.injected_signal_fraction, 100 * c.relative_gain, marker="o")
            ax.axhline(0, linewidth=1)
            ax.set_xlabel("Injected fine-signal fraction")
            ax.set_ylabel("Recovered partial-pooling gain (%)")
            ax.set_title("Signal detectability calibration")
            fig.tight_layout()
            p = outdir / "signal_calibration.png"
            fig.savefig(p, dpi=180, bbox_inches="tight")
            plt.close(fig)
            plots.append(p)

        if extra.get("permutation") is not None:
            pm = extra["permutation"]
            vals = pm.null_table.partial_pooling_gain.to_numpy()
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.hist(vals, bins=30)
            ax.axvline(pm.observed_partial_pooling_gain, linewidth=2)
            ax.set_xlabel("Partial-pooling gain")
            ax.set_ylabel("Permutation count")
            ax.set_title("Common-taxonomy permutation null")
            fig.tight_layout()
            p = outdir / "permutation_null.png"
            fig.savefig(p, dpi=180, bbox_inches="tight")
            plt.close(fig)
            plots.append(p)

        return plots

    def _write_html(self, outdir, fit, extra, plots, manifest):
        def table_html(df, selected=False):
            return dataframe_to_html_table(
                df,
                selected_level=(fit.selected_level if selected else None),
                level_column="level",
            )

        sections = [
            "<h1>TransferRes analysis report</h1>",
            f"<p><b>Software:</b> TransferRes {APP_VERSION}</p>",
            f"<p><b>Selected transferable resolution:</b> {fit.selected_level}</p>",
            "<h2>Resolution table</h2>",
            table_html(fit.level_table, selected=True),
            "<h2>Transition table</h2>",
            table_html(fit.transition_table),
        ]

        if extra.get("bootstrap") is not None:
            sections += [
                "<h2>Biological-unit bootstrap</h2>",
                table_html(extra["bootstrap"].raw_level_table),
                table_html(extra["bootstrap"].selection_frequency),
            ]
        if extra.get("ident") is not None:
            sections += [
                "<h2>Structural identifiability</h2>",
                table_html(extra["ident"].summary_table),
                table_html(extra["ident"].support_frontier),
            ]
        if extra.get("metric") is not None:
            sections += ["<h2>Metric sensitivity</h2>", table_html(extra["metric"].interpretation_table)]
        if extra.get("feature_sensitivity") is not None:
            sections += ["<h2>Feature-count sensitivity</h2>", table_html(extra["feature_sensitivity"].interpretation_table)]
        if extra.get("permutation") is not None:
            p = extra["permutation"]
            sections += [
                "<h2>Common-taxonomy permutation</h2>",
                f"<p>Partial gain P(greater) = {p.p_partial_greater:.6g}; "
                f"two-sided P = {p.p_partial_two_sided:.6g}</p>",
            ]
        if extra.get("calibration") is not None:
            sections += [
                "<h2>Signal-injection calibration</h2>",
                table_html(extra["calibration"].calibration_table),
            ]

        sections += ["<h2>Software-assisted interpretation</h2>", "<ul>" + "".join(f"<li>{x}</li>" for x in build_interpretation(fit, extra)) + "</ul>"]

        sections.append("<h2>Plots</h2>")
        for p in plots:
            sections.append(
                f'<div><img src="{p.name}" style="max-width:900px;width:100%;"></div>'
            )

        sections += [
            "<h2>Interpretation contract</h2>",
            "<ul>"
            "<li>Resolution selection is defined by minimum raw held-out SSE.</li>"
            "<li>Partial-pooling gains are diagnostics, not the discrete selection rule.</li>"
            "<li>Zero transferable gain does not imply absence of finer biological structure.</li>"
            "<li>Bootstrap and perturbation-bootstrap intervals must be labeled by their resampling unit.</li>"
            "</ul>",
        ]

        style = """
        :root{
          --border:#d8dbe2;
          --header:#f3f5f8;
          --selected:#eef6ff;
          --text:#20242a;
          --muted:#68707b;
        }
        body{
          font-family:Segoe UI,Arial,sans-serif;
          max-width:1120px;
          margin:34px auto;
          padding:0 22px 60px;
          color:var(--text);
          line-height:1.45;
        }
        h1{font-size:30px;margin-bottom:8px}
        h2{font-size:21px;margin-top:34px;margin-bottom:10px}
        p{margin:7px 0}
        .table-scroll{
          width:100%;
          overflow-x:auto;
          margin:12px 0 28px 0;
          border:1px solid var(--border);
          border-radius:8px;
          background:white;
        }
        table.tr-table{
          border-collapse:separate;
          border-spacing:0;
          width:max-content;
          min-width:100%;
          font-size:13px;
        }
        table.tr-table th,
        table.tr-table td{
          padding:8px 10px;
          border-right:1px solid var(--border);
          border-bottom:1px solid var(--border);
          white-space:nowrap;
          text-align:right;
          vertical-align:top;
        }
        table.tr-table th{
          background:var(--header);
          font-weight:600;
          position:sticky;
          top:0;
          z-index:1;
        }
        table.tr-table th:first-child,
        table.tr-table td:first-child{
          text-align:left;
        }
        table.tr-table th:last-child,
        table.tr-table td:last-child{
          border-right:none;
        }
        table.tr-table tbody tr:last-child td{
          border-bottom:none;
        }
        table.tr-table tbody tr.selected-row td{
          background:var(--selected);
          font-weight:600;
        }
        img{
          display:block;
          max-width:100%;
          height:auto;
          margin:12px 0 30px 0;
          border:1px solid #e4e6eb;
          border-radius:8px;
        }
        ul{padding-left:22px}
        """

        doc = (
            "<!doctype html><html><head><meta charset='utf-8'>"
            f"<style>{style}</style></head><body>{''.join(sections)}</body></html>"
        )
        (outdir / "report.html").write_text(doc, encoding="utf-8")

    # ---------------- results ----------------
    def _populate_results(self, fit, outdir, extra, plots):
        self.result_selected.set(fit.selected_level)
        self.result_path.set(str(outdir))
        self.result_text.delete("1.0", "end")
        self.interpretation_box.delete("1.0","end")
        for line in build_interpretation(fit, extra):
            self.interpretation_box.insert("end", "• "+line+"\n")
        self.result_text.insert("end", "RESOLUTION TABLE\n")
        self.result_text.insert("end", fit.level_table.to_string(index=False))
        self.result_text.insert("end", "\n\nTRANSITION TABLE\n")
        self.result_text.insert("end", fit.transition_table.to_string(index=False))

        if extra.get("metric") is not None:
            self.result_text.insert("end", "\n\nMETRIC INTERPRETATION\n")
            self.result_text.insert(
                "end", extra["metric"].interpretation_table.to_string(index=False)
            )
        if extra.get("feature_sensitivity") is not None:
            self.result_text.insert("end", "\n\nFEATURE-COUNT INTERPRETATION\n")
            self.result_text.insert(
                "end",
                extra["feature_sensitivity"].interpretation_table.to_string(index=False),
            )
        if extra.get("permutation") is not None:
            p = extra["permutation"]
            self.result_text.insert(
                "end",
                f"\n\nPERMUTATION\nPartial one-sided P = {p.p_partial_greater:.6g}\n"
                f"Partial two-sided P = {p.p_partial_two_sided:.6g}\n",
            )

        self.plot_list.delete(0, "end")
        for p in plots:
            self.plot_list.insert("end", p.name)

    # ---------------- open actions ----------------

    def _show_about(self):
        messagebox.showinfo(
            "About TransferRes",
            "TransferRes 0.6.0-rc1\n\n"
            "Transferable perturbational resolution analysis.\n\n"
            "Primary selection rule: minimum raw held-out SSE.\n"
            "Nested feature selection is recommended.\n"
            "Partial pooling is diagnostic, not the discrete level-selection rule."
        )

    def _open_output(self):
        p = Path(self.output_dir.get())
        p.mkdir(parents=True, exist_ok=True)
        try:
            os.startfile(str(p))
        except Exception:
            messagebox.showinfo("Output folder", str(p))

    def _open_last_result(self):
        if not self.last_output:
            messagebox.showinfo("Results", "No analysis has been run.")
            return
        try:
            os.startfile(str(self.last_output))
        except Exception:
            messagebox.showinfo("Results folder", str(self.last_output))

    def _open_html_report(self):
        if not self.last_output:
            messagebox.showinfo("Results", "No analysis has been run.")
            return
        p = self.last_output / "report.html"
        if p.exists():
            webbrowser.open(p.as_uri())
        else:
            messagebox.showinfo("Report", "No HTML report found.")

    def _open_error_log(self):
        if not self.last_output:
            messagebox.showinfo("Error log", "No analysis output exists.")
            return
        p = self.last_output / "error_log.txt"
        if not p.exists():
            messagebox.showinfo("Error log", "No error log exists for the last run.")
            return
        try:
            os.startfile(str(p))
        except Exception:
            webbrowser.open(p.as_uri())

    def _open_selected_plot(self):
        if not self.last_output:
            return
        sel = self.plot_list.curselection()
        if not sel:
            return
        p = self.last_output / self.plot_list.get(sel[0])
        try:
            os.startfile(str(p))
        except Exception:
            webbrowser.open(p.as_uri())

    def _log(self, text):
        self.status.insert("end", text)
        self.status.see("end")
        self.update_idletasks()


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(run_executable_self_test())
    TransferResApp().mainloop()
