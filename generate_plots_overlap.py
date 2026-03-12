import scripts.plots as plots
# import scripts.utils as utils

# ==========================================
#               CONSTANTS
# ==========================================
# TSV col name to human readable name
overlap_mapping = {
    "overlap_seq_abs_gx_attn": "|gx| & attn",
    "overlap_seq_abs_ig_attn": "|ig| & attn",
    "overlap_gx_pos": "gx > 0",
    "overlap_gx_abs": "|gx|",
    "overlap_gx_abs_diag": "|gx (Diagonal)|",
    "overlap_gx_abs_anti": "|gx (Anti-Diag.)|",
    "overlap_ig_abs": "|ig|",
    "overlap_ig_pos": "ig > 0",
    "overlap_ig_abs_diag": "|ig (Diagonal)|",
    "overlap_ig_abs_anti": "|ig (Anti-Diag.)|",
}

# Generate plot for specific GO terms 
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/overlap_summary_specific.tsv",
    output_dir="filtered_metrics_r1/top_10_overlap/specific",
    metric_mapping=overlap_mapping,
    y_limit=120 
)

#  Generate plot for all GO terms including propagated
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/overlap_summary_propagated.tsv",
    output_dir="filtered_metrics_r1/top_10_overlap/propagated",
    metric_mapping=overlap_mapping,
    y_limit=120
)