import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

META = "fold_switching_metadata.tsv"
PROP_DIR = "raw_data/models_propagated_data"
OUT_PNG = "plots/delta_vs_rmsd.png"
# OUT_PNG = OUT_DIR / "delta_vs_TM.png"
# OUT_PNG = OUT_DIR / "delta_vs_seqid.png"

MODELS = [("fusion", "pred_prob", 0.1), ("struct", "struct_prob", 0.4), ("esm", "esm_prob", 0.1)]
MODEL_COLS = {
    "pred_prob": ("pred_prob_fold1", "pred_prob_fold2", "fusion_delta"),
    "struct_prob": ("struct_prob_fold1", "struct_prob_fold2", "struct_delta"),
    "esm_prob": ("esm_prob_fold1", "esm_prob_fold2", "esm_delta"),
}
COLS = ["go_id", "pred_prob", "struct_prob", "esm_prob"]
TM_COL = "TM-score"
ID_COL = "identity"
RMSD_COL = "rmsd"
# markers cycle if more proteins than entries
MARKERS = "os^vDPhP*X"


def filter_m_by_model_thresholds(m: pd.DataFrame) -> pd.DataFrame:
    """Keep GO rows where ≥1 model has fold1 or fold2 ≥ thr; NaN that model's cols otherwise."""
    any_active = pd.Series(False, index=m.index)
    for _, col, thr in MODELS:
        c1, c2 = f"{col}_fold1", f"{col}_fold2"
        active = (m[c1] >= thr) | (m[c2] >= thr)
        any_active |= active
        k1, k2, kd = MODEL_COLS[col]
        m.loc[~active, [k1, k2, kd]] = float("nan")
    return m.loc[any_active].copy()


meta = pd.read_csv(META, sep="\t")


records: list[dict] = []
for _, row in meta.iterrows():
    protein, f1, f2 = row["Protein"], str(row["fold1"]), str(row["fold2"])
    tm = float(row[TM_COL])
    id = float(row[ID_COL])
    rmsd = float(row[RMSD_COL])
    p1 = pd.read_csv(PROP_DIR / f"{f1}.csv", usecols=COLS)
    p2 = pd.read_csv(PROP_DIR / f"{f2}.csv", usecols=COLS)
    m = p1.merge(p2, on="go_id", how="outer", suffixes=("_fold1", "_fold2"))
    m["fusion_delta"] = m["pred_prob_fold2"] - m["pred_prob_fold1"]
    m["struct_delta"] = m["struct_prob_fold2"] - m["struct_prob_fold1"]
    m["esm_delta"] = m["esm_prob_fold2"] - m["esm_prob_fold1"]
    m = filter_m_by_model_thresholds(m)
    if m.empty:
        continue
    for label, col, _ in MODELS:
        kd = MODEL_COLS[col][2]
        for v in m[kd].abs().dropna():
            # records.append({"Protein": protein, "model": label, "abs_delta": float(v), "tm": tm})
            # records.append({"Protein": protein, "model": label, "abs_delta": float(v), "id": id})
            records.append({"Protein": protein, "model": label, "abs_delta": float(v), "rmsd": rmsd})

if not records:
    raise SystemExit("No points to plot (empty after filters).")

long_df = pd.DataFrame.from_records(records)
proteins = list(dict.fromkeys(meta["Protein"]))  # stable order as in meta
prot_marker = {p: MARKERS[i % len(MARKERS)] for i, p in enumerate(proteins)}
model_color = {label: f"C{i}" for i, (label, _, _) in enumerate(MODELS)}

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
for (protein, model), g in long_df.groupby(["Protein", "model"], sort=False):
    ax.scatter(
        # g["tm"],
        # g["id"],
        g["rmsd"],
        g["abs_delta"],
        s=14,
        alpha=1,
        marker=prot_marker[protein],
        c=model_color[model],
        edgecolors="none",
    )

ax.set_ylabel("|Δ Probability| (Fold2 − Fold1)")
# ax.set_xlabel(TM_COL)
ax.set_xlabel("RMSD")
# ax.set_xlabel("Sequence identity")
# ax.set_title("|Δ Probability| vs. pair TM-score")
ax.set_title("|Δ Probability| vs. pair RMSD")
# ax.set_title("|Δ Probability| vs. pair sequence identity")
ax.grid(True, alpha=0.25)
ax.set_ylim(bottom=0)
# ax.set_xlim(0.85, 1.05) # For seq identity
# ax.set_xlim(0.2, 0.6) # For TM score
ax.set_xlim(2.2, 4) # For RMSD

model_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=model_color[m], markersize=9, label=m, linestyle="none")
    for m, _, _ in MODELS
]
prot_handles = [
    Line2D([0], [0], marker=prot_marker[p], color="0.25", markerfacecolor="0.25", markersize=8, label=p, linestyle="none")
    for p in proteins
    if p in long_df["Protein"].unique()
]
leg_m = ax.legend(handles=model_handles, title="model", loc="upper left", bbox_to_anchor=(1, 0.5), frameon=True)
ax.add_artist(leg_m)
ax.legend(handles=prot_handles, title="Protein", loc="lower left", bbox_to_anchor=(1, 0.5), frameon=True)

fig.savefig(OUT_PNG, dpi=150)
plt.close(fig)
