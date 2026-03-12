import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score
import scripts.utils as utils

def run_auroc_evaluation(biolip_path, output_tsv_path, specific_terms_only):
    
    # --- 1. Load and Parse Ground Truth (BioLiP) ---
    if not os.path.exists(biolip_path):
        print(f"Error: Ground truth file not found at {biolip_path}. Skipping...")
        return
        
    df_biolip = pd.read_csv(biolip_path, sep='\t')
    
    ground_truth = {}
    for _, row in df_biolip.iterrows():
        uniprot_id = str(row['uniprot_id']).strip()
        ligand_id = str(row['ligand_id']).strip()
        go_terms_raw = str(row['go_terms_pruned'])
        res_str = str(row['binding_resnums_uniprot_mapped'])
        
        if pd.isna(go_terms_raw) or go_terms_raw == 'nan' or go_terms_raw == 'N/A':
            continue
        if pd.isna(res_str) or res_str == 'nan' or not res_str.strip():
            continue
            
        res_list = [int(x) for x in res_str.split(',')]
        
        for gt in go_terms_raw.split(','):
            go_term = gt.strip().replace(':', '_')
            key = (uniprot_id, ligand_id, go_term)
            if key not in ground_truth:
                ground_truth[key] = set()
            ground_truth[key].update(res_list)

    # --- 2. Directory processing setup ---
    directory = "local"
    
    # Use utils to fetch the requested files
    seq_files = utils.get_target_files(directory, specific_terms_only=specific_terms_only)
    
    results = []
    print(f"Found {len(seq_files)} sequence files to evaluate.")

    for seq_file in seq_files:
        struct_file = seq_file.replace("local_seq_", "local_struct_")
        
        if not os.path.exists(struct_file):
            continue
            
        # Use utils to extract metadata
        protein, csv_go_term = utils.get_protein_and_go(seq_file)
        uniprot_id = protein.split('-')[1]
        
        matching_keys = [k for k in ground_truth.keys() if k[0] == uniprot_id and k[2] == csv_go_term]
        if not matching_keys:
            continue

        # Use utils to load, clean, align, and add derived columns
        df_seq, df_struct, fair_max_limit = utils.load_and_prepare_data(seq_file, struct_file)
        
        # Ensure dataframes are sorted by residue_idx to prevent random array shuffling 
        df_seq = df_seq.sort_values('residue_idx')
        df_struct = df_struct.sort_values('residue_idx')

        for key in matching_keys:
            _, ligand_id, _ = key
            
            # Get valid true residues within truncated sequence length
            raw_true_residues = ground_truth[key]
            true_residues = {res for res in raw_true_residues if res <= fair_max_limit}
            
            if len(true_residues) == 0:
                continue
                
            # Robust Helper function to calculate AUROC
            def get_auroc(df, score_col):
                if score_col not in df.columns:
                    return None
                    
                # Dynamically build y_true to match EXACTLY the rows present in this dataframe
                y_true = np.array([1 if r in true_residues else 0 for r in df['residue_idx']])
                
                # AUROC is undefined if there are no positive examples OR no negative examples
                if len(np.unique(y_true)) < 2:
                    return None
                    
                y_pred = df[score_col].values
                
                # Catch edge case where the model output an array of all exactly equal values
                if len(np.unique(y_pred)) == 1:
                    return 0.5 
                    
                return roc_auc_score(y_true, y_pred)

            # Calculate and append all metrics (Sequence and Structure)
            results.append({
                "Protein": protein,
                "UniProt": uniprot_id,
                "Ligand": ligand_id,
                "GO_Term": csv_go_term,
                "Valid_Sequence_Length": fair_max_limit, 
                "Ground_Truth_Count": len(true_residues),
                
                # Sequence Metrics (6)
                "Seq_signed_AUROC": get_auroc(df_seq, 'gradxinput'),
                "Seq_pos_AUROC": get_auroc(df_seq, 'pos_gradxinput'),
                "Seq_abs_AUROC": get_auroc(df_seq, 'abs_gradxinput'),
                "Seq_ig_signed_AUROC": get_auroc(df_seq, 'ig_signed'),
                "Seq_ig_pos_AUROC": get_auroc(df_seq, 'ig_pos'),
                "Seq_ig_abs_AUROC": get_auroc(df_seq, 'abs_ig'),
                
                # Structure Metrics (14)
                "Struct_signed_AUROC": get_auroc(df_struct, 'gx_signed'),
                "Struct_pos_AUROC": get_auroc(df_struct, 'gx_pos'),
                "Struct_abs_AUROC": get_auroc(df_struct, 'gx_abs'),
                "Struct_abs_anti_AUROC": get_auroc(df_struct, 'gx_abs_anti'),
                "Struct_abs_diag_AUROC": get_auroc(df_struct, 'gx_abs_diag'),
                "Struct_gx_signed_diag_AUROC": get_auroc(df_struct, 'gx_signed_diag'),
                "Struct_gx_signed_anti_AUROC": get_auroc(df_struct, 'gx_signed_anti'),
                "Struct_ig_signed_AUROC": get_auroc(df_struct, 'ig_signed'),
                "Struct_ig_pos_AUROC": get_auroc(df_struct, 'ig_pos'),
                "Struct_ig_abs_AUROC": get_auroc(df_struct, 'ig_abs'),
                "Struct_ig_signed_diag_AUROC": get_auroc(df_struct, 'ig_signed_diag'),
                "Struct_ig_signed_anti_AUROC": get_auroc(df_struct, 'ig_signed_anti'),
                "Struct_ig_abs_diag_AUROC": get_auroc(df_struct, 'ig_abs_diag'),
                "Struct_ig_abs_anti_AUROC": get_auroc(df_struct, 'ig_abs_anti')
            })

    # --- 3. Save to TSV ---
    if results:
        output_df = pd.DataFrame(results)
        os.makedirs(os.path.dirname(output_tsv_path) or ".", exist_ok=True)
        output_df.to_csv(output_tsv_path, sep='\t', index=False)
        print(f"Processed {len(results)} matches. Results saved to: {output_tsv_path}")
    else:
        print("No matching ground truth records were found.")


# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Run for Specific Terms ONLY
run_auroc_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_specificGO_combined_residues_ligand_GOterm.tsv",
    output_tsv_path="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_specificGOterm_AUROC.tsv",
    specific_terms_only=True
)

# 2. Run for ALL Terms (Propagated)
run_auroc_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_propagatedGO_combined_residues_ligand_GOterm.tsv", 
    output_tsv_path="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_propagatedGOterm_AUROC.tsv",
    specific_terms_only=False
)