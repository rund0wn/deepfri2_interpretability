import scripts.plots as plots

fi_mapping = {
    # "Seq_signed_FI": "Seq: grad x input",
    "Seq_pos_FI": "Seq: (grad x input) > 0",
    "Seq_abs_FI": "Seq: |grad x input|",
    # "Seq_ig_signed_FI": "Seq: ig",
    "Seq_ig_pos_FI": "Seq: ig > 0",
    "Seq_ig_abs_FI": "Seq: |ig|",

    # "Struct_signed_FI": "grad x input",
    "Struct_pos_FI": "(grad x input) > 0",
    "Struct_abs_FI": "|grad x input|",
    "Struct_abs_anti_FI": "|grad x input| (Anti-Diag.)",
    "Struct_abs_diag_FI": "|grad x input| (Diagonal)",
    # "Struct_gx_signed_diag_FI": "grad x input (Diagonal)",
    # "Struct_gx_signed_anti_FI": "grad x input (Anti-Diag.)",
    # "Struct_ig_signed_FI": "ig",
    "Struct_ig_pos_FI": "ig > 0",
    "Struct_ig_abs_FI": "|ig|",
    "Struct_ig_abs_anti_FI": "|ig| (Anti-Diag.)",
    "Struct_ig_abs_diag_FI": "|ig| (Diagonal)",
    # "Struct_ig_signed_diag_FI": "ig (Diagonal)",
    # "Struct_ig_signed_anti_FI": "ig (Anti-Diag.)"
}

plots.generate_fi_per_protein_plots(
    input_file="integrated_gradients/auprc_fi_specific.tsv",
    output_dir="integrated_gradients/auprc_fi_plots/per_protein_fi/specific",
    metric_mapping=fi_mapping
)

plots.generate_fi_per_protein_plots(
    input_file="integrated_gradients/auprc_fi_propagated.tsv",
    output_dir="integrated_gradients/auprc_fi_plots/per_protein_fi/propagated",
    metric_mapping=fi_mapping
)