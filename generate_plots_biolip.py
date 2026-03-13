import scripts.plots as plots

# 1. Define the exact columns you want to plot and their legend labels
biolip_mapping = {
    # "Seq_signed_overlap": "seq: grad x input",
    "DeepFRI_1.0_sal_overlap": "DeepFRI 1: saliency",
    "Seq_attn_overlap": "seq: attn",

    "Seq_pos_overlap": "seq: (gx) > 0",
    "Seq_abs_overlap": "seq: |gx|",
    # "Seq_ig_signed_overlap": "seq: ig",
    "Seq_ig_pos_overlap": "seq: ig > 0",
    "Seq_ig_abs_overlap": "seq: |ig|",
    # "Struct_signed_overlap": "grad x input",
    "Struct_pos_overlap": "(gx) > 0",
    "Struct_abs_overlap": "|gx|",
    "Struct_abs_anti_overlap": "|gx| (Anti-Diag.)",
    "Struct_abs_diag_overlap": "|gx| (Diagonal)",
    # "Struct_signed_diag_overlap": "grad x input (Diagonal)",
    # "Struct_signed_anti_overlap": "grad x input (Anti-Diag.)",
    # "Struct_ig_signed_overlap": "ig",
    "Struct_ig_pos_overlap": "ig > 0",
    "Struct_ig_abs_overlap": "|ig|",
    # "Struct_ig_signed_diag_overlap": "ig (Diagonal)",
    # "Struct_ig_signed_anti_overlap": "ig (Anti-Diag.)",
    "Struct_ig_abs_diag_overlap": "|ig| (Diagonal)",
    "Struct_ig_abs_anti_overlap": "|ig| (Anti-Diag.)"
}


# 2. Run for specific GO terms
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/biolip_multimetric_overlap_combined_residues_ligand_specificGOterm.tsv",
    output_dir="filtered_metrics_r1/top_10_biolip/specific",
    metric_mapping=biolip_mapping,
    y_limit=120
)

# 3. Run for propagated GO terms
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/biolip_multimetric_overlap_combined_residues_ligand_propagatedGOterm.tsv",
    output_dir="filtered_metrics_r1/top_10_biolip/propagated",
    metric_mapping=biolip_mapping,
    y_limit=120
)
