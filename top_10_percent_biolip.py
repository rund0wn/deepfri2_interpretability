import pandas as pd
import os
import scripts.utils as utils
import numpy as np

def run_biolip_evaluation(biolip_path, output_tsv_path, specific_terms_only):
    
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
            go_term_biolip = gt.strip().replace(':', '_')
            key = (uniprot_id, ligand_id, go_term_biolip)
            if key not in ground_truth:
                ground_truth[key] = set()
            ground_truth[key].update(res_list)

    # --- 2. Directory processing setup ---
    directory = "local"
    all_files = utils.get_target_files(directory, specific_terms_only=specific_terms_only, deepfri=True)
    seq_files = [f for f in all_files if "local_seq_" in f]

    results = []
    print(f"Found {len(all_files)} files to evaluate.")

    
    for seq_file in seq_files:
        struct_file = seq_file.replace("local_seq_", "local_struct_")
        deepfri_file = seq_file.replace("local_seq_", "deepfri_v1_").replace("v4", "v6") # Making my life easier, but they are actually v6
        
        # Get protein and GO term pair, then check if there is BioLip ground truth for them
        protein, go_term = utils.get_protein_and_go(seq_file)
        matching_keys = [k for k in ground_truth.keys() if k[0] == protein and k[2] == go_term]
        if not matching_keys:
            continue
        
        print(f'Protein: {protein}, GO Term: {go_term}')
        if os.path.exists(deepfri_file):
            df_seq, df_struct, df_deepfri, fair_max_limit = utils.load_and_prepare_data(seq_file, struct_file, deepfri_file)
        else:
            df_seq, df_struct, fair_max_limit = utils.load_and_prepare_data(seq_file, struct_file)
            df_deepfri = None

        # Calculate top 10% based on the shared, identical length
        top_k = max(1, int(fair_max_limit * 0.10))

        for key in matching_keys:
            _, ligand_id, _ = key
            raw_true_residues = ground_truth[key]
            true_residues = {res for res in raw_true_residues if res <= fair_max_limit}
            n_true = len(true_residues)
            if n_true == 0:
                continue


            def get_gt_overlap(df, col):
                # if df is None or col not in df.columns:
                if df is None:
                    return np.nan
                top_residues = set(df.nlargest(top_k, col)['residue_idx'])
                return (len(top_residues.intersection(true_residues)) / n_true) * 100
            
            results.append({
                "Protein": protein,
                # "UniProt": uniprot_id,
                "Ligand": ligand_id,
                "GO_Term": go_term,
                "Valid_Sequence_Length": fair_max_limit, 
                "Ground_Truth_Count": n_true,

                # Added DeepFRI 1.0 
                "DeepFRI_1.0_sal_overlap": get_gt_overlap(df_deepfri, 'saliency'),
                
                # # Sequence Metrics
                "Seq_attn_overlap": get_gt_overlap(df_seq, 'attn'),
                # "Seq_signed_overlap": get_gt_overlap(df_seq, 'gradxinput'),
                "Seq_pos_overlap": get_gt_overlap(df_seq, 'pos_gradxinput'),
                "Seq_abs_overlap": get_gt_overlap(df_seq, 'abs_gradxinput'),
                # "Seq_ig_signed_overlap": get_gt_overlap(df_seq, 'ig_signed'),
                "Seq_ig_pos_overlap": get_gt_overlap(df_seq, 'ig_pos'),
                "Seq_ig_abs_overlap": get_gt_overlap(df_seq, 'abs_ig'),
                
                # # Structure Metrics
                # "Struct_signed_overlap": get_gt_overlap(df_struct, 'gx_signed'),
                "Struct_pos_overlap": get_gt_overlap(df_struct, 'gx_pos'),
                "Struct_abs_overlap": get_gt_overlap(df_struct, 'gx_abs'),
                "Struct_abs_anti_overlap": get_gt_overlap(df_struct, 'gx_abs_anti'),
                "Struct_abs_diag_overlap": get_gt_overlap(df_struct, 'gx_abs_diag'),
                # "Struct_signed_diag_overlap": get_gt_overlap(df_struct, 'gx_signed_diag'),
                # "Struct_signed_anti_overlap": get_gt_overlap(df_struct, 'gx_signed_anti'),
                # "Struct_ig_signed_overlap": get_gt_overlap(df_struct, 'ig_signed'),
                "Struct_ig_pos_overlap": get_gt_overlap(df_struct, 'ig_pos'),
                "Struct_ig_abs_overlap": get_gt_overlap(df_struct, 'ig_abs'),
                # "Struct_ig_signed_diag_overlap": get_gt_overlap(df_struct, 'ig_signed_diag'),
                # "Struct_ig_signed_anti_overlap": get_gt_overlap(df_struct, 'ig_signed_anti'),
                "Struct_ig_abs_diag_overlap": get_gt_overlap(df_struct, 'ig_abs_diag'),
                "Struct_ig_abs_anti_overlap": get_gt_overlap(df_struct, 'ig_abs_anti')
            })

    # 3. Save to TSV
    if results:
        output_df = pd.DataFrame(results)
        
        # Create the output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_tsv_path), exist_ok=True)
        
        output_df.to_csv(output_tsv_path, sep='\t', index=False, na_rep='NaN')
        print(f"Processed {len(results)} matches. Results saved to: {output_tsv_path}")
    else:
        print("No matching ground truth records were found.")
            

# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Run for Specific Terms ONLY
run_biolip_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_specificGO_combined_residues_ligand_GOterm.tsv",
    output_tsv_path="filtered_metrics_r1/biolip_multimetric_overlap_combined_residues_ligand_specificGOterm.tsv",
    specific_terms_only=True
)

# 2. Run for ALL Terms (Propagated)
run_biolip_evaluation(
    biolip_path="biolip_mappings/biolip_eval_subset_propagatedGO_combined_residues_ligand_GOterm.tsv",
    output_tsv_path="filtered_metrics_r1/biolip_multimetric_overlap_combined_residues_ligand_propagatedGOterm.tsv",
    specific_terms_only=False
)