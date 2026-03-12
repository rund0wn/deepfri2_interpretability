import scripts.plots as plots
# import scripts.utils as utils

# ==========================================
#               CONSTANTS
# ==========================================
OVERLAP_COMPARISONS = [
        # ("overlap_gx_signed",      "gradxinput",     "gx_signed"),
        ("overlap_seq_abs_gx_attn", "attn",           "abs_gradxinput"),
        ("overlap_seq_abs_ig_attn", "attn",           "abs_ig"),

        ("overlap_seq_abs_ig_attn", "attn",           "abs_ig"),
        ("overlap_gx_pos",         "pos_gradxinput", "gx_pos"),
        ("overlap_gx_abs",         "abs_gradxinput", "gx_abs"),
        # ("overlap_gx_signed_diag", "gradxinput",     "gx_signed_diag"),
        # ("overlap_gx_signed_anti", "gradxinput",     "gx_signed_anti"),
        ("overlap_gx_abs_diag",    "abs_gradxinput", "gx_abs_diag"),
        ("overlap_gx_abs_anti",    "abs_gradxinput", "gx_abs_anti"),
        # ("overlap_ig_signed",      "ig_signed",      "ig_signed"),
        ("overlap_ig_abs",         "abs_ig",         "ig_abs"),
        ("overlap_ig_pos",         "ig_pos",         "ig_pos"),
        # ("overlap_ig_signed_diag", "ig_signed",      "ig_signed_diag"),
        # ("overlap_ig_signed_anti", "ig_signed",      "ig_signed_anti"),
        ("overlap_ig_abs_diag",    "abs_ig",         "ig_abs_diag"), 
        ("overlap_ig_abs_anti",    "abs_ig",         "ig_abs_anti")
    ]

def format_label(col_name):
    """Converts raw column names like 'gx_abs_diag' to '|gx_abs| (Diagonal)'."""
    base = col_name.replace("_diag", "").replace("_anti", "")
    label = f"|{base}|" if "abs" in base else base
    
    if "_diag" in col_name:
        label += " (Diagonal)"
    elif "_anti" in col_name:
        label += " (Anti-Diag.)"
        
    return label

# 1. Build the dynamic mapping dictionary using utils
overlap_mapping = {
    col_name: format_label(struct_col) 
    # for col_name, seq_col, struct_col in utils.OVERLAP_COMPARISONS
    for col_name, seq_col, struct_col, in OVERLAP_COMPARISONS
}

# 2. Generate plot for specific GO terms 
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/overlap_summary_specific.tsv",
    output_dir="filtered_metrics_r1/top_10_overlap/specific",
    metric_mapping=overlap_mapping,
    y_limit=120 
)

# 3. Generate plot for all GO terms including propagated
plots.generate_grouped_barplot(
    input_file="filtered_metrics_r1/overlap_summary_propagated.tsv",
    output_dir="filtered_metrics_r1/top_10_overlap/propagated",
    metric_mapping=overlap_mapping,
    y_limit=120
)