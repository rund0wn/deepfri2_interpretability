import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import textwrap

def get_go_dict(seq_path, struct_path):
    """Extracts GO ID -> GO Name mapping from prediction files."""
    go_dict = {}
    for path in [seq_path, struct_path]:
        if os.path.exists(path):
            tmp = pd.read_csv(path)
            # The prediction files use "GO:", which matches your standardized df
            go_dict.update(dict(zip(tmp["GO_term"], tmp["GO_name"])))
    return go_dict

def format_go_label(go_id, go_dict):
    """Combines GO ID and wrapped GO name for the X-axis label."""
    name = go_dict.get(go_id, "Unknown")
    wrapped_name = "\n".join(textwrap.wrap(name, width=22))
    return f"{go_id}\n{wrapped_name}"

def generate_grouped_barplot(input_file, output_dir, metric_mapping, y_limit=120):
    df_ss = pd.read_csv(input_file, sep="\t")
    df_ss["GO_Term"] = df_ss["GO_Term"].str.replace("GO_", "GO:")
    # df_ss["Protein"] = df_ss["Protein"].str.strip()

    seq_pred_path = "/Users/rund/deepfri2_local/interpretability_local/prediction_scores/seq_model_predictions.csv"
    struct_pred_path = "/Users/rund/deepfri2_local/interpretability_local/prediction_scores/struct_model_predictions.csv"

    # 1. Build the dictionary of names
    go_dict = get_go_dict(seq_pred_path, struct_pred_path)

    if os.path.exists(seq_pred_path) and os.path.exists(struct_pred_path):
        df_seq_pred = pd.read_csv(seq_pred_path)
        df_struct_pred = pd.read_csv(struct_pred_path)

        df_seq_pred = df_seq_pred.rename(columns={"GO_term": "GO_Term", "predicted_p>=0.30": "seq_P"})
        df_struct_pred = df_struct_pred.rename(columns={"GO_term": "GO_Term", "predicted_p>=0.30": "struct_P"})
        
        df_seq_pred["Protein"] = df_seq_pred["Protein"].str.split("-").str[1]
        df_struct_pred["Protein"] = df_struct_pred["Protein"].str.split("-").str[1]

        df_ss = df_ss.merge(df_seq_pred[["Protein", "GO_Term", "seq_P"]], on=["Protein", "GO_Term"], how="left")
        df_ss = df_ss.merge(df_struct_pred[["Protein", "GO_Term", "struct_P"]], on=["Protein", "GO_Term"], how="left")
    else:
        df_ss["seq_P"] = None
        df_ss["struct_P"] = None


    # 2. Apply the mapping to create a new label column
    df_ss["GO_Label"] = df_ss["GO_Term"].apply(lambda x: format_go_label(x, go_dict))

    # df_ss["UniProt"] = df_ss["Protein"].str.split("-").str[1]
    id_vars = ["Protein", "GO_Label", "seq_P", "struct_P"]

    plot_metrics = [col for col in metric_mapping.keys() if col in df_ss.columns]

    df_melted_ss = df_ss.melt(
        id_vars=id_vars,
        value_vars=plot_metrics,
        var_name="Metric",
        value_name="Overlap (%)",
    )

    df_melted_ss["Metric"] = df_melted_ss["Metric"].map(metric_mapping)

    # (Keep existing color palette code here...)
    seq_metrics = [col for col in metric_mapping.keys() if col.startswith('Seq_') and col in plot_metrics]
    struct_metrics = [col for col in metric_mapping.keys() if col.startswith('Struct_') and col in plot_metrics]

    # Check if we have specialized categories
    if seq_metrics or struct_metrics:
        # Build the dual-color palette (Blues and Oranges)
        seq_display_names = [metric_mapping[col] for col in seq_metrics]
        struct_display_names = [metric_mapping[col] for col in struct_metrics]
        
        seq_colors = sns.color_palette("Blues", n_colors=max(3, len(seq_display_names)))[:len(seq_display_names)]
        struct_colors = sns.color_palette("Oranges", n_colors=max(3, len(struct_display_names)))[:len(struct_display_names)]
        
        palette = {}
        for name, color in zip(seq_display_names, seq_colors): palette[name] = color
        for name, color in zip(struct_display_names, struct_colors): palette[name] = color
    else:
        # Fallback to viridis if no "Seq_" or "Struct_" prefixes are found
        palette = "viridis"

    os.makedirs(output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")

    for protein, group in df_melted_ss.groupby("Protein"):
        n_terms = group["GO_Label"].nunique()
        fig_width = max(8, n_terms * 2.5)
        fig, ax = plt.subplots(figsize=(fig_width, 8)) # Taller to fit names

        # Plot using GO_Label
        sns.barplot(data=group, x="GO_Label", y="Overlap (%)", hue="Metric", palette=palette, ax=ax, errorbar=None)

        # Annotations using GO_Label
        go_terms = group["GO_Label"].unique()
        for i, term in enumerate(go_terms):
            term_data = group[group["GO_Label"] == term].iloc[0]
            s_val, st_val = term_data["seq_P"], term_data["struct_P"]

            s_txt = f"{s_val:.2f}" if pd.notnull(s_val) else "N/A"
            st_txt = f"{st_val:.2f}" if pd.notnull(st_val) else "N/A"

            ax.text(
                i, y_limit * 0.875, f"Seq: {s_txt}\nStr: {st_txt}",
                ha="center", va="bottom", fontsize=9, fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1),
            )

        ax.set_ylim(0, y_limit)
        ax.set_title(f"{protein}", fontsize=14, fontweight="bold", pad=20)
        ax.set_xlabel("") 
        
        # Style the X tick labels
        plt.xticks(fontsize=9, fontweight="medium")
        ax.set_ylabel("Overlap (%)")

        # Legend
        ax.legend(
            loc="lower center", bbox_to_anchor=(0.5, 1.12),
            ncol=min(5, len(group["Metric"].unique())),
            fontsize=10, title="Metric", title_fontsize=11, frameon=True,
        )

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{protein}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"Plots saved to {output_dir}")

def generate_fi_per_protein_plots(input_file, output_dir, metric_mapping):
    """
    Generates per-protein grouped barplots specifically for AUPRC Fold-Improvement.
    Draws a red dashed line at y=1.0 to indicate random baseline.
    Includes seq and struct prediction scores as annotations.
    """
    df = pd.read_csv(input_file, sep="\t")
    
    # Load prediction scores (optional - skip if files don't exist)
    seq_pred_path = "/Users/rund/deepfri2_local/interpretability_local/prediction_scores/seq_model_predictions.csv"
    struct_pred_path = "/Users/rund/deepfri2_local/interpretability_local/prediction_scores/struct_model_predictions.csv"

    if os.path.exists(seq_pred_path) and os.path.exists(struct_pred_path):
        df_seq_pred = pd.read_csv(seq_pred_path)
        df_struct_pred = pd.read_csv(struct_pred_path)
        df_seq_pred["GO_term"] = df_seq_pred["GO_term"].str.replace("GO:", "GO_")
        df_struct_pred["GO_term"] = df_struct_pred["GO_term"].str.replace("GO:", "GO_")
  
        df_seq_pred = df_seq_pred.rename(
            columns={"GO_term": "GO_Term", "predicted_p>=0.30": "seq_P"}
        )
        df_seq_pred["Protein"] = df_seq_pred["Protein"].str.strip()
        
        df_struct_pred = df_struct_pred.rename(
            columns={"GO_term": "GO_Term", "predicted_p>=0.30": "struct_P"}
        )
        df_struct_pred["Protein"] = df_struct_pred["Protein"].str.strip()

        df = df.merge(
            df_seq_pred[["Protein", "GO_Term", "seq_P"]],
            on=["Protein", "GO_Term"],
            how="left",
        )

        df = df.merge(
            df_struct_pred[["Protein", "GO_Term", "struct_P"]],
            on=["Protein", "GO_Term"],
            how="left",
        )
    else:
        # Create columns with None if files don't exist
        df["seq_P"] = None
        df["struct_P"] = None
    
    # Filter only columns that exist in mapping
    plot_metrics = [col for col in metric_mapping.keys() if col in df.columns]
    
    # Melt the dataframe - preserve seq_P and struct_P through melting
    df_melted = df.melt(
        id_vars=["Protein", "UniProt", "GO_Term", "seq_P", "struct_P"],
        value_vars=plot_metrics,
        var_name="Metric",
        value_name="Fold Improvement"
    )
    
    df_melted["Metric"] = df_melted["Metric"].map(metric_mapping)
    
    os.makedirs(output_dir, exist_ok=True)
    sns.set_theme(style="whitegrid")
    
    # Group by protein to make one plot per protein
    for protein, group in df_melted.groupby("Protein"):
        n_terms = group["GO_Term"].nunique()
        fig_width = max(8, n_terms * 2.5) 
        
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        
        sns.barplot(
            data=group,
            x="GO_Term",
            y="Fold Improvement",
            hue="Metric",
            palette="viridis",
            ax=ax,
            errorbar=None
        )
        
        # Draw the random baseline at FI = 1.0
        ax.axhline(y=1.0, color='red', linestyle='--', linewidth=2, zorder=10, label='Random Baseline (1.0x)')
        
        # Dynamic Y-axis: Add a little headroom above the highest bar
        max_fi = group["Fold Improvement"].max()
        y_limit = max(2.0, max_fi * 1.15)
        ax.set_ylim(0, y_limit)
        
        # Add Text Annotations
        go_terms = group["GO_Term"].unique()
        for i, term in enumerate(go_terms):
            term_data = group[group["GO_Term"] == term].iloc[0]
            s_val, st_val = term_data["seq_P"], term_data["struct_P"]

            s_txt = f"{s_val:.2f}" if pd.notnull(s_val) else "N/A"
            st_txt = f"{st_val:.2f}" if pd.notnull(st_val) else "N/A"

            # Position text relative to the dynamic y_limit
            text_y_pos = y_limit * 0.875

            ax.text(
                i,
                text_y_pos,
                f"Seq: {s_txt}\nStr: {st_txt}",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1),
            )
        
        ax.set_title(f"{protein} - AUPRC Fold Improvement", fontsize=14, fontweight="bold", pad=20)
        ax.set_xlabel("GO Term", fontsize=12, fontweight="bold")
        ax.set_ylabel("Fold Improvement (over random)", fontsize=12, fontweight="bold")
        
        # Legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=min(4, len(handles)), frameon=False)
        
        plt.tight_layout()
        outfile = os.path.join(output_dir, f"{protein}_FI.png")
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"Per-protein FI plots saved to {output_dir}")