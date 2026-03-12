import pandas as pd
import glob
import os

# ==========================================
# 1. CONSTANTS
# ==========================================
# List of specific terms for current protein set
TERM_LIST = [
    'GO_0000030', 'GO_0046872', 'GO_0005254', 'GO_0051087', 'GO_0042802', 
    'GO_0000981', 'GO_0003677', 'GO_0005515', 'GO_0008270', 'GO_0019899', 
    'GO_0050897', 'GO_0005524', 'GO_0031625', 'GO_0035251', 'GO_0020037', 
    'GO_0005506'
]

OVERLAP_COMPARISONS = [
        # ("overlap_gx_signed",      "gradxinput",     "gx_signed"),
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

# ==========================================
# 2. METADATA EXTRACTION
# ==========================================
def get_protein_and_go(filepath):
    """Extracts protein and GO term from your specific filename format."""
    basename = os.path.basename(filepath.split("/")[-1])

    # Works for seq, struct, and deepfri_v1
    if basename.startswith("local_"):
        name_parts = basename.replace("local_seq_", "").replace("local_struct_", "").replace(".csv", "").split("__")
        protein = name_parts[0].split("-")[1]
        go_term = name_parts[1]
    else:
        name_parts = basename.replace("deepfri_v1_", "").replace(".csv", "").split("__")
        protein = name_parts[0].split("-")[1]
        go_term = name_parts[1]
    
    return protein, go_term

# ==========================================
# 3. FILE FILTERING
# ==========================================
def get_target_files(directory, specific_terms_only=True, deepfri=False):
    """Finds seq files and optionally filters them by TERM_LIST."""
    seq_files = glob.glob(os.path.join(directory, "local_seq_*.csv"))
    struct_files = glob.glob(os.path.join(directory, "local_struct_*.csv"))
    all_files = seq_files + struct_files

    if deepfri:
        deepfri_files = glob.glob(os.path.join(directory, "deepfri_v1_*.csv"))
        all_files = all_files + deepfri_files

    if not specific_terms_only:
        return all_files
        
    filtered_files = []
    for file in all_files:
        _, go_term = get_protein_and_go(file)
        if go_term in TERM_LIST:
            filtered_files.append(file)
    
    return filtered_files

# ==========================================
# 4. CORE DATA PROCESSING
# ==========================================
def load_and_prepare_data(seq_file, struct_file, deepfri_file=None):
    """Loads, cleans, truncates, and adds derived columns for a file pair."""
    # 1. Load the data
    df_seq = pd.read_csv(seq_file)
    df_struct = pd.read_csv(struct_file)

    # Make fair comparison (truncate to remove trailing residues)
    df_seq = df_seq.loc[(df_seq.drop('residue_idx', axis=1) != 0).any(axis=1)].copy()
    df_struct = df_struct.loc[(df_struct.drop('residue_idx', axis=1) != 0).any(axis=1)].copy()

    # set residue_idx for struct to start at 1
    df_struct.index = df_struct.index + 1
    df_struct['residue_idx'] = df_struct.index

    # Calculate limits
    max_seq_idx = df_seq['residue_idx'].max()
    max_struct_idx = df_struct['residue_idx'].max()
    fair_max_limit = min(max_seq_idx, max_struct_idx)


    # --- OPTIONAL DEEPFRI PROCESSING ---
    if deepfri_file:
        df_deepfri = pd.read_csv(deepfri_file)
        df_deepfri = df_deepfri.loc[(df_deepfri.drop('residue_idx', axis=1) != 0).any(axis=1)].copy()
        df_deepfri.index = df_deepfri.index + 1
        
        max_deepfri_idx = df_deepfri['residue_idx'].max()
        fair_max_limit = min(fair_max_limit, max_deepfri_idx)


    if deepfri_file:
        df_deepfri = df_deepfri[df_deepfri['residue_idx'] <= fair_max_limit]

    # Trim to the fair_max_limit
    df_struct = df_struct[df_struct['residue_idx'] <= fair_max_limit]
    df_seq = df_seq[df_seq['residue_idx'] <= fair_max_limit]

    # 2. Add the positive 'gradxinput' column to the dataframes
    df_seq['pos_gradxinput'] = df_seq['gradxinput'].clip(lower=0)
    df_struct['gx_pos'] = df_struct['gx_signed'].clip(lower=0)
    df_seq['ig_pos'] = df_seq['ig_signed'].clip(lower=0)
    df_struct['ig_pos'] = df_struct['ig_signed'].clip(lower=0)

    # 3. Add the absolute 'gradxinput' column to the sequence dataframe
    df_seq['abs_gradxinput'] = df_seq['gradxinput'].abs()
    df_seq['abs_ig'] = df_seq['ig_signed'].abs()

    if deepfri_file:
        return df_seq, df_struct, df_deepfri, fair_max_limit
    else:
        return df_seq, df_struct, fair_max_limit


# ### TEST
# directory = "../local"

# # Set specific_terms_only to False if you want all files
# target_files = get_target_files(directory,  specific_terms_only=False, deepfri=True)
# seq_files = [f for f in target_files if "local_seq_" in f]

# for seq_file in seq_files:
#     struct_file = seq_file.replace("local_seq_", "local_struct_")
#     deepfri_file = seq_file.replace("local_seq_", "deepfri_v1_").replace("v4", "v6") # Making my life easier, but they are actually v6

#     protein, go_term = get_protein_and_go(seq_file)

#     if os.path.exists(deepfri_file):
#         df_seq, df_struct, df_deepfri, fair_max_limit = load_and_prepare_data(seq_file, struct_file, deepfri_file)
#     else:
#         df_seq, df_struct, fair_max_limit = load_and_prepare_data(seq_file, struct_file)


