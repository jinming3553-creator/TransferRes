
from pathlib import Path
import sys
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

from TransferRes_GUI import dataframe_to_html_table

df=pd.DataFrame({
    "level":["broad","class_name"],
    "heldout_sse":[110.539,69.3552],
    "gain_vs_coarsest":[0.0,0.372573],
    "raw_incremental_gain_vs_parent":[None,0.372573],
    "n_biological_units":[8,8],
})

h=dataframe_to_html_table(
    df,
    selected_level="class_name",
    level_column="level",
)

for col in df.columns:
    token=f"<th>{col}</th>"
    assert token in h, f"Missing independent header cell: {col}"

assert h.count("<th>")==len(df.columns)
assert h.count("</th>")==len(df.columns)
assert 'class="selected-row"' in h
assert 'class="table-scroll"' in h
assert '<td>class_name</td>' in h

out=ROOT/"_html_smoke_output.html"
out.write_text(
    "<html><body>"+h+"</body></html>",
    encoding="utf-8",
)
print("HTML_RENDERER_SMOKE_TEST=PASS")
print(f"N_HEADERS={len(df.columns)}")
print(out)
