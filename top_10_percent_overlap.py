import pandas as pd
import os
import sys
import scripts.utils as utils

directory = "local"

def process_files(seq_files, output_name):
    results = []
    print(f"Processing {len(seq_files)} files for {output_name}...")

    for seq_file in seq_files:
        # File path
        struct_file = seq_file.replace("local_seq_", "local_struct_")

        if os.path.exists(struct_file):
            # get protein/GO term, and prepare dfs
            protein, go_term = utils.get_protein_and_go(seq_file)
            df_seq, df_struct, fair_max_limit = utils.load_and_prepare_data(seq_file, struct_file)

            # Determine how many residues make up 10%
            top_k = max(1, int(fair_max_limit * 0.10))

            # Helper function: compute overlap between seq and struct dfs
            def get_overlap(seq_col, struct_col):
                seq_top = set(df_seq.nlargest(top_k, seq_col)["residue_idx"])
                struct_top = set(df_struct.nlargest(top_k, struct_col)["residue_idx"])
                return (len(seq_top.intersection(struct_top)) / top_k) * 100
            
            # Helper function: compute overlap within a single df (added for seq attn)
            def get_internal_overlap(df, col1, col2):
                if col1 in df.columns and col2 in df.columns:
                    top1 = set(df.nlargest(top_k, col1)["residue_idx"])
                    top2 = set(df.nlargest(top_k, col2)["residue_idx"])
                    return (len(top1.intersection(top2)) / top_k) * 100
                return None

            # Initialize result dict with labels
            res_dict = {
                "Protein": protein,
                "GO_Term": go_term
            }

            # Get overlap of top 10% between attn and other seq metrics (added for seq)
            res_dict["overlap_seq_abs_gx_attn"] = get_internal_overlap(df_seq, "attn", "abs_gradxinput")
            res_dict["overlap_seq_abs_ig_attn"] = get_internal_overlap(df_seq, "attn", "abs_ig")

            # Dynamically calculate overlaps based on the mapping
            for col_name, seq_col, struct_col in utils.OVERLAP_COMPARISONS:
                # Ensure the columns actually exist before trying to calculate
                if seq_col in df_seq.columns and struct_col in df_struct.columns:
                    res_dict[col_name] = get_overlap(seq_col, struct_col)
                else:
                    res_dict[col_name] = None


            results.append(res_dict)


    if results:
        output_df = pd.DataFrame(results)
        output_df.to_csv(output_name, sep="\t", index=False)
        print(f"Processed {len(results)} pairs. Results saved to {output_name}.")
    else:
        print(f"No results found for {output_name}.\n")


# # ==========================================
# # MAIN EXECUTION
# # ==========================================

# # 1. Run for specific terms ONLY
specific_files = utils.get_target_files(directory, specific_terms_only=True)
seq_specific_files = [f for f in specific_files if "local_seq_" in f]
process_files(seq_specific_files, "filtered_metrics_r1/overlap_summary_specific.tsv")

# # 2. Run for ALL terms
all_files = utils.get_target_files(directory, specific_terms_only=False)
seq_all_files = [f for f in all_files if "local_seq_" in f]
process_files(seq_all_files, "filtered_metrics_r1/overlap_summary_propagated.tsv")