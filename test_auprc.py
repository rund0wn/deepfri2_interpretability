import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

# import pandas as pd
# import numpy as np
# import glob
# import os
# from sklearn.metrics import average_precision_score

# # --- 1. Load and Parse Ground Truth (BioLiP) ---
# biolip_path = "biolip_eval_subset_specificGO_combined_residues_ligand_GOterm.tsv" 
# df_biolip = pd.read_csv(biolip_path, sep='\t')

# ground_truth = {}
# for _, row in df_biolip.iterrows():
#     uniprot_id = str(row['uniprot_id']).strip()
#     ligand_id = str(row['ligand_id']).strip()
#     go_terms_raw = str(row['go_terms_pruned'])
#     res_str = str(row['binding_resnums_uniprot_mapped'])
    
#     if pd.isna(go_terms_raw) or go_terms_raw == 'nan' or go_terms_raw == 'N/A':
#         continue
#     if pd.isna(res_str) or res_str == 'nan' or not res_str.strip():
#         continue
        
#     res_list = [int(x) for x in res_str.split(',')]
    
#     for gt in go_terms_raw.split(','):
#         go_term = gt.strip().replace(':', '_')
#         key = (uniprot_id, ligand_id, go_term)
#         if key not in ground_truth:
#             ground_truth[key] = set()
#         ground_truth[key].update(res_list)

# # --- 2. Directory processing setup ---
# directory = "/mnt/vdb2/var/storage/deepfri2/training/notebooks/interpretability/local/"
# seq_files = glob.glob(os.path.join(directory, "local_seq_*.csv"))

# results = []

# for seq_file in seq_files:
#     struct_file = seq_file.replace("local_seq_", "local_struct_")
    
#     if not os.path.exists(struct_file):
#         continue
        
#     basename = os.path.basename(seq_file)
#     name_parts = basename.replace("local_seq_", "").replace(".csv", "").split("__")
#     protein = name_parts[0]
#     csv_go_term = name_parts[1]
#     uniprot_id = protein.split('-')[1]
    
#     matching_keys = [k for k in ground_truth.keys() if k[0] == uniprot_id and k[2] == csv_go_term]
#     if not matching_keys:
#         continue

#     # Load data
#     df_seq_raw = pd.read_csv(seq_file)
#     df_struct_raw = pd.read_csv(struct_file)
    
#     # --- STRICT FAIRNESS TRUNCATION LOGIC ---
#     seq_features = ['attn', 'cos', 'gradxinput']
#     struct_features = ['gx_signed', 'gx_abs', 'gx_signed_diag', 'gx_signed_anti', 'gx_abs_diag', 'gx_abs_anti']
    
#     df_seq_valid = df_seq_raw.loc[(df_seq_raw[seq_features] != 0.0).any(axis=1)].copy()
#     df_struct_valid = df_struct_raw.loc[(df_struct_raw[struct_features] != 0.0).any(axis=1)].copy()
    
#     if df_seq_valid.empty or df_struct_valid.empty:
#         continue
        
#     max_seq_idx = df_seq_valid['residue_idx'].max()
#     max_struct_idx = df_struct_valid['residue_idx'].max() + 1
    
#     fair_max_limit = min(max_seq_idx, max_struct_idx)
    
#     df_seq = df_seq_valid[df_seq_valid['residue_idx'] <= fair_max_limit].copy()
#     df_struct = df_struct_valid[df_struct_valid['residue_idx'] < fair_max_limit].copy()
    
#     df_seq['abs_gradxinput'] = df_seq['gradxinput'].abs()
    
#     # Extract prediction arrays (ensure sorted by residue_idx so they map to y_true correctly)
#     df_seq = df_seq.sort_values('residue_idx')
#     df_struct = df_struct.sort_values('residue_idx')
    
#     # ----------------------------------------

#     for key in matching_keys:
#         _, ligand_id, _ = key
        
#         # Build the y_true binary "Answer Key" array
#         y_true = np.zeros(fair_max_limit, dtype=int)
        
#         raw_true_residues = ground_truth[key]
#         for res in raw_true_residues:
#             # Check bounds and convert 1-based res to 0-based array index
#             if 1 <= res <= fair_max_limit:
#                 y_true[res - 1] = 1 
                
#         # If there are no positive examples within bounds, AUPRC is undefined
#         n_true = sum(y_true)
#         if n_true == 0:
#             continue
            
#         # The Baseline is the score of the "blindfolded toddler" (random guessing)
#         # It's simply the fraction of true binding residues vs total residues
#         baseline_auprc = n_true / fair_max_limit
            
#         # Helper function to calculate AUPRC
#         def get_auprc(df, score_col):
#             y_pred = df[score_col].values
            
#             # Catch edge case where the model output an array of all exactly equal values
#             if len(np.unique(y_pred)) == 1:
#                 return baseline_auprc
                
#             return average_precision_score(y_true, y_pred)

#         # Calculate AUPRC metrics
#         auprc_seq_signed = get_auprc(df_seq, 'gradxinput')
#         auprc_seq_abs = get_auprc(df_seq, 'abs_gradxinput')
        
#         auprc_struct_signed = get_auprc(df_struct, 'gx_signed')
#         auprc_struct_abs = get_auprc(df_struct, 'gx_abs')
#         auprc_struct_abs_anti = get_auprc(df_struct, 'gx_abs_anti')
#         auprc_struct_abs_diag = get_auprc(df_struct, 'gx_abs_diag')
        
#         results.append({
#             "Protein": protein,
#             "UniProt": uniprot_id,
#             "Ligand": ligand_id,
#             "GO_Term": csv_go_term,
#             "Valid_Sequence_Length": fair_max_limit, 
#             "Ground_Truth_Count": n_true,
#             "Random_Baseline_AUPRC": baseline_auprc, # Added so you can see what score to beat
#             "Seq_signed_AUPRC": auprc_seq_signed,
#             "Seq_abs_AUPRC": auprc_seq_abs,
#             "Struct_signed_AUPRC": auprc_struct_signed,
#             "Struct_abs_AUPRC": auprc_struct_abs,
#             "Struct_abs_anti_AUPRC": auprc_struct_abs_anti,
#             "Struct_abs_diag_AUPRC": auprc_struct_abs_diag
#         })

# # --- 3. Save to TSV ---
# output_df = pd.DataFrame(results)
# output_tsv_path = "biolip_multimetric_overlap_combined_residues_ligand_GOterm_AUPRC.tsv"
# output_df.to_csv(output_tsv_path, sep='\t', index=False)
# print(f"Results saved to: {output_tsv_path}")

##### FIGURE
# --- TARGET TO INSPECT ---
target_protein = "AF-P00452-F1-model_v4_A"
target_ligand = "A1A3L"
target_go_term = "GO_0005515"
# -------------------------

# --- 1. Load and Parse Ground Truth (BioLiP) ---
biolip_path = "biolip_eval_subset_specificGO_combined_residues_ligand_GOterm.tsv" 
df_biolip = pd.read_csv(biolip_path, sep='\t')

ground_truth = {}
for _, row in df_biolip.iterrows():
    uniprot_id = str(row['uniprot_id']).strip()
    ligand_id = str(row['ligand_id']).strip()
    go_terms_raw = str(row['go_terms_pruned'])
    res_str = str(row['binding_resnums_uniprot_mapped'])
    
    if pd.isna(go_terms_raw) or go_terms_raw == 'nan' or pd.isna(res_str):
        continue
        
    res_list = [int(x) for x in res_str.split(',')]
    for gt in go_terms_raw.split(','):
        key = (uniprot_id, ligand_id, gt.strip().replace(':', '_'))
        if key not in ground_truth:
            ground_truth[key] = set()
        ground_truth[key].update(res_list)

# --- 2. Process the Specific Target ---
directory = "/mnt/vdb2/var/storage/deepfri2/training/notebooks/interpretability/local/"
seq_file = os.path.join(directory, f"local_seq_{target_protein}__{target_go_term}.csv")
struct_file = os.path.join(directory, f"local_struct_{target_protein}__{target_go_term}.csv")

uniprot_id = target_protein.split('-')[1]
key = (uniprot_id, target_ligand, target_go_term)

if not os.path.exists(seq_file) or not os.path.exists(struct_file) or key not in ground_truth:
    print("Error: Files or ground truth not found for this target.")
    exit()

df_seq_raw = pd.read_csv(seq_file)
df_struct_raw = pd.read_csv(struct_file)

# Fairness Truncation
seq_features = ['attn', 'cos', 'gradxinput']
struct_features = ['gx_signed', 'gx_abs', 'gx_signed_diag', 'gx_signed_anti', 'gx_abs_diag', 'gx_abs_anti']

df_seq_valid = df_seq_raw.loc[(df_seq_raw[seq_features] != 0.0).any(axis=1)].copy()
df_struct_valid = df_struct_raw.loc[(df_struct_raw[struct_features] != 0.0).any(axis=1)].copy()

fair_max_limit = min(df_seq_valid['residue_idx'].max(), df_struct_valid['residue_idx'].max() + 1)

df_seq = df_seq_valid[df_seq_valid['residue_idx'] <= fair_max_limit].copy().sort_values('residue_idx')
df_struct = df_struct_valid[df_struct_valid['residue_idx'] < fair_max_limit].copy().sort_values('residue_idx')
df_seq['abs_gradxinput'] = df_seq['gradxinput'].abs()

# Build y_true
y_true = np.zeros(fair_max_limit, dtype=int)
for res in ground_truth[key]:
    if 1 <= res <= fair_max_limit:
        y_true[res - 1] = 1 

# --- 3. Calculate and Plot PR Curves ---
plt.figure(figsize=(10, 8))

# Define the metrics we want to plot
metrics_to_plot = {
    'Sequence (Absolute)': (df_seq['abs_gradxinput'].values, 'purple'),
    'Structure (Absolute)': (df_struct['gx_abs'].values, 'teal'),
    'Structure (Anti-Diag)': (df_struct['gx_abs_anti'].values, 'orange')
}

for label, (y_pred, color) in metrics_to_plot.items():
    # Calculate Precision, Recall, and Average Precision (AUPRC)
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    auprc = average_precision_score(y_true, y_pred)
    
    # Plot the curve
    plt.plot(recall, precision, color=color, lw=2, label=f'{label} (AUPRC = {auprc:.2f})')

# Baseline (Random Guessing Performance)
baseline = sum(y_true) / len(y_true)
plt.plot([0, 1], [baseline, baseline], linestyle='--', color='gray', label=f'Random Baseline ({baseline:.3f})')

plt.xlabel('Recall (Fraction of True Binding Sites Found)', fontsize=12)
plt.ylabel('Precision (Fraction of Predictions That Are Correct)', fontsize=12)
plt.title(f'Precision-Recall Curve: {target_protein} ({target_ligand})', fontsize=14, pad=15)
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, alpha=0.3)

output_image = "auprc_comparison_test.png"
plt.savefig(output_image, dpi=300, bbox_inches='tight')
print(f"AUPRC Curve successfully saved to {output_image}")
