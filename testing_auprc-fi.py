import pandas as pd
import numpy as np
import os
from sklearn.metrics import average_precision_score
import scripts.utils as utils

def run_auprc_fi_evaluation(biolip_path, output_tsv_path, specific_terms_only):
    
    # --- 1. Load and Parse Ground Truth (BioLiP) ---
    if not os.path.exists(biolip_path):
        print(f"Error: File {biolip_path} not found.")
        return
        
    df_biolip = pd.read_csv(biolip_path, sep='\t')
    ground_truth = {}
    
    for _, row in df_biolip.iterrows():
        uniprot_id = str(row['uniprot_id']).strip()
        ligand_id = str(row['ligand_id']).strip()
        go_terms_raw = str(row['go_terms_pruned'])
        res_str = str(row['binding_resnums_uniprot_mapped'])
        
        if pd.isna(go_terms_raw) or go_terms_raw == 'nan' or pd.isna(res_str) or res_str == 'nan':
            continue
            
        res_list = [int(x) for x in res_str.split(',')]
        
        for gt in go_terms_raw.split(','):
            go_term = gt.strip().replace(':', '_')
            key = (uniprot_id, go_term)
            if key not in ground_truth:
                ground_truth[key] = set()
            ground_truth[key].update(res_list)

    # --- 2. Directory processing setup ---
    directory = "local"
    seq_files = utils.get_target_files(directory, specific_terms_only=specific_terms_only)
    results = []

    for seq_file in seq_files:
        struct_file = seq_file.replace("local_seq_", "local_struct_")
        if not os.path.exists(struct_file): continue
            
        protein, csv_go_term = utils.get_protein_and_go(seq_file)
        uniprot_id = protein.split('-')[1]
        
        key = (uniprot_id, csv_go_term)
        if key not in ground_truth: continue

        # Use utils to load, clean, and align
        df_seq, df_struct, fair_max_limit = utils.load_and_prepare_data(seq_file, struct_file)
        
        # Sort values to ensure y_true and y_pred perfectly align
        df_seq = df_seq.sort_values('residue_idx')
        df_struct = df_struct.sort_values('residue_idx')

        # Build y_true
        y_true = np.zeros(fair_max_limit, dtype=int)
        for res in ground_truth[key]:
            if 1 <= res <= fair_max_limit:
                y_true[res - 1] = 1 
                
        n_true = sum(y_true)
        if n_true == 0: continue
            
        baseline_auprc = n_true / fair_max_limit
        
        def get_auprc_fi(df, score_col):
            if score_col not in df.columns: return None
            y_pred = df[score_col].values
            if len(np.unique(y_pred)) == 1:
                return 1.0 # Random baseline (FI = 1.0)
            auprc = average_precision_score(y_true, y_pred)
            return auprc / baseline_auprc 

        # Calculate and append all metrics (Sequence and Structure)
        results.append({
                "Protein": protein,
                "UniProt": uniprot_id,
                "Ligand": ligand_id,
                "GO_Term": csv_go_term,
                "Valid_Sequence_Length": fair_max_limit, 
                "Ground_Truth_Count": n_true,
                
                # Sequence Metrics (6)
                "Seq_signed_FI": get_auprc_fi(df_seq, 'gradxinput'),
                "Seq_pos_FI": get_auprc_fi(df_seq, 'pos_gradxinput'),
                "Seq_abs_FI": get_auprc_fi(df_seq, 'abs_gradxinput'),
                "Seq_ig_signed_FI": get_auprc_fi(df_seq, 'ig_signed'),
                "Seq_ig_pos_FI": get_auprc_fi(df_seq, 'ig_pos'),
                "Seq_ig_abs_FI": get_auprc_fi(df_seq, 'abs_ig'),
                
                # Structure Metrics (14)
                "Struct_signed_FI": get_auprc_fi(df_struct, 'gx_signed'),
                "Struct_pos_FI": get_auprc_fi(df_struct, 'gx_pos'),
                "Struct_abs_FI": get_auprc_fi(df_struct, 'gx_abs'),
                "Struct_abs_anti_FI": get_auprc_fi(df_struct, 'gx_abs_anti'),
                "Struct_abs_diag_FI": get_auprc_fi(df_struct, 'gx_abs_diag'),
                "Struct_gx_signed_diag_FI": get_auprc_fi(df_struct, 'gx_signed_diag'),
                "Struct_gx_signed_anti_FI": get_auprc_fi(df_struct, 'gx_signed_anti'),
                "Struct_ig_signed_FI": get_auprc_fi(df_struct, 'ig_signed'),
                "Struct_ig_pos_FI": get_auprc_fi(df_struct, 'ig_pos'),
                "Struct_ig_abs_FI": get_auprc_fi(df_struct, 'ig_abs'),
                "Struct_ig_signed_diag_FI": get_auprc_fi(df_struct, 'ig_signed_diag'),
                "Struct_ig_signed_anti_FI": get_auprc_fi(df_struct, 'ig_signed_anti'),
                "Struct_ig_abs_diag_FI": get_auprc_fi(df_struct, 'ig_abs_diag'),
                "Struct_ig_abs_anti_FI": get_auprc_fi(df_struct, 'ig_abs_anti')
            })

    # --- 3. Save to TSV ---
    if results:
        output_df = pd.DataFrame(results)
        os.makedirs(os.path.dirname(output_tsv_path) or ".", exist_ok=True)
        output_df.to_csv(output_tsv_path, sep='\t', index=False)
        print(f"Processed {len(results)} matches. Saved to: {output_tsv_path}")

# Run for specific terms
run_auprc_fi_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_specificGO_combined_residues_ligand_GOterm.tsv",
    output_tsv_path="integrated_gradients/auprc_fi_specific.tsv",
    specific_terms_only=True
)

# Run for propagated terms
run_auprc_fi_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_propagatedGO_combined_residues_ligand_GOterm.tsv",
    output_tsv_path="integrated_gradients/auprc_fi_propagated.tsv",
    specific_terms_only=False
)