import pandas as pd
from utils.information_content.calc_MF_IC_sem_sim import calculate_IC

## Input files/paths
META = "fold_switching_metadata.tsv"
PROP_DIR = "raw_data/models_propagated_data"
DEEPFRI1_DIR = "raw_data/deepfriv1_data"
OUT_DELTA = "results/results_fold_switching_deltas.tsv"
OUT_CROSS = "results/results_fold_switching_threshold_crossings.tsv"
OUT_JACCARD = "results/results_fold_switching_jaccard.tsv"

# (model label, csv column, threshold)
MODELS = [("fusion", "pred_prob", 0.1), ("struct", "struct_prob", 0.4), ("esm", "esm_prob", 0.1), ("deepfri1", "deepfri1_prob", 0.29)]
# score column prefix on merged frame -> (fold1, fold2, delta) output cols for delta TSV
DELTA_COLS = [
    ("pred_prob", 0.1, ("fusion_fold1", "fusion_fold2", "fusion_delta")),
    ("struct_prob", 0.4, ("struct_fold1", "struct_fold2", "struct_delta")),
    ("esm_prob", 0.1, ("esm_fold1", "esm_fold2", "esm_delta")),
    ("deepfri1_prob", 0.29, ("deepfri1_fold1", "deepfri1_fold2", "deepfri1_delta")),
]
COLS = ["go_id", "pred_prob", "struct_prob", "esm_prob"]


## Calculate Jaccard similarity between two sets of GO terms
def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def main() -> None:
    meta = pd.read_csv(META, sep="\t")

    delta_parts: list[pd.DataFrame] = []
    cross_rows: list[dict] = []
    jacc_rows: list[dict] = []

    for _, row in meta.iterrows():
        protein, f1, f2 = row["Protein"], str(row["fold1"]), str(row["fold2"])

        # Models
        p1 = pd.read_csv(PROP_DIR + f"/{f1}.csv", usecols=COLS)
        p2 = pd.read_csv(PROP_DIR + f"/{f2}.csv", usecols=COLS)
        m = p1.merge(p2, on="go_id", how="outer", suffixes=("_fold1", "_fold2"))

        # DeepFRI v1
        d1 = pd.read_csv(DEEPFRI1_DIR + f"/{f1}_MF_predictions.csv", skiprows = 1, usecols = ['GO_term/EC_number', 'Score'])
        d1.rename(columns={'GO_term/EC_number': 'go_id', 'Score': 'deepfri1_prob'}, inplace = True)
        d2 = pd.read_csv(DEEPFRI1_DIR + f"/{f2}_MF_predictions.csv", skiprows = 1, usecols = ['GO_term/EC_number', 'Score'])
        d2.rename(columns={'GO_term/EC_number': 'go_id', 'Score': 'deepfri1_prob'}, inplace = True)
        m2 = d1.merge(d2, on="go_id", how="outer", suffixes=("_fold1", "_fold2"))

        m.insert(0, "Protein", protein)
        m.insert(1, "fold1", f1)
        m.insert(2, "fold2", f2)
        # m["go_term_name"] = m["go_id"].map(lambda g: names.get(g, ""))

        m["fusion_fold1"] = m["pred_prob_fold1"]
        m["fusion_fold2"] = m["pred_prob_fold2"]
        m["fusion_delta"] = m["fusion_fold2"] - m["fusion_fold1"]
        m["struct_fold1"] = m["struct_prob_fold1"]
        m["struct_fold2"] = m["struct_prob_fold2"]
        m["struct_delta"] = m["struct_fold2"] - m["struct_fold1"]
        m["esm_fold1"] = m["esm_prob_fold1"]
        m["esm_fold2"] = m["esm_prob_fold2"]
        m["esm_delta"] = m["esm_fold2"] - m["esm_fold1"]

        # Add DeepFRI v1 scores and deltas to m
        m["deepfri1_prob_fold1"] = m2["deepfri1_prob_fold1"]
        m["deepfri1_prob_fold2"] = m2["deepfri1_prob_fold2"]
        m["deepfri1_fold1"] = m2["deepfri1_prob_fold1"]
        m["deepfri1_fold2"] = m2["deepfri1_prob_fold2"]
        m["deepfri1_delta"] = m["deepfri1_fold2"] - m["deepfri1_fold1"]

        # Add Information Content (IC) of GO terms
        m["go_ic"] = m["go_id"].map(calculate_IC)

        # Delta TSV: only report per-model scores+delta if that model is >= threshold on fold1 or fold2
        d = m[
            [
                "Protein",
                "fold1",
                "fold2",
                "go_id",
                "go_ic",
                # "go_term_name",
                "fusion_fold1",
                "fusion_fold2",
                "fusion_delta",
                "struct_fold1",
                "struct_fold2",
                "struct_delta",
                "esm_fold1",
                "esm_fold2",
                "esm_delta",
                "deepfri1_fold1",
                "deepfri1_fold2",
                "deepfri1_delta",
            ]
        ].copy()
        any_active = pd.Series(False, index=m.index)
        for col, thr, (k1, k2, kd) in DELTA_COLS:
            c1, c2 = f"{col}_fold1", f"{col}_fold2"
            active = (m[c1] >= thr) | (m[c2] >= thr)
            any_active |= active
            d.loc[~active, [k1, k2, kd]] = float("nan")
        delta_parts.append(d.loc[any_active])

        for model, col, thr in MODELS:
            c1, c2 = f"{col}_fold1", f"{col}_fold2"
            s1, s2 = m[c1], m[c2]
            pass1, pass2 = s1 >= thr, s2 >= thr
            set1 = set(m.loc[pass1, "go_id"])
            set2 = set(m.loc[pass2, "go_id"])
            inter, o1, o2 = set1 & set2, set1 - set2, set2 - set1
            jacc_rows.append(
                {
                    "Protein": protein,
                    "fold1": f1,
                    "fold2": f2,
                    "model": model,
                    # "n_fold1": len(set1),
                    # "n_fold2": len(set2),
                    # "n_shared": len(inter),
                    # "n_only_fold1": len(o1),
                    # "n_only_fold2": len(o2),
                    "jaccard": jaccard(set1, set2),
                }
            )

            xmask = pass1 != pass2
            for _, r in m.loc[xmask].iterrows():
                p1, p2 = bool(r[c1] >= thr), bool(r[c2] >= thr)
                direction = "below_to_above" if not p1 and p2 else "above_to_below"
                cross_rows.append(
                    {
                        "Protein": protein,
                        "fold1": f1,
                        "fold2": f2,
                        "model": model,
                        "go_id": r["go_id"],
                        # "go_term_name": r["go_term_name"],
                        "score_fold1": r[c1],
                        "score_fold2": r[c2],
                        "delta": r[c2] - r[c1],
                        "pass_fold1": p1,
                        "pass_fold2": p2,
                        "direction": direction,
                    }
                )

    d_all = pd.concat(delta_parts, ignore_index=True)
    d_all.to_csv(OUT_DELTA, sep="\t", index=False)
    pd.DataFrame(cross_rows).to_csv(OUT_CROSS, sep="\t", index=False)
    pd.DataFrame(jacc_rows).to_csv(OUT_JACCARD, sep="\t", index=False)
    n_cross = len(cross_rows)
    print(f"Wrote {OUT_DELTA} ({len(d_all)} rows), {OUT_CROSS} ({n_cross} rows), {OUT_JACCARD}.")


if __name__ == "__main__":
    main()
