import pandas as pd

# Paths
deepfri_dir = '../local'
saliency_json = deepfri_dir + '/deepfri_v1/sanity_MF_saliency_maps_v2.json'

# Read the saliency maps JSON file
saliency_df = pd.read_json(saliency_json, orient='index')

# Create a separate file for each protein-GO term pair saliency values
def create_seperate_files(df, output_dir):
    for protein, data in saliency_df.iterrows():
        GO_ids = data['GO_ids']
        saliency_maps = data['saliency_maps']
        id_saliency_map = dict(zip(GO_ids, saliency_maps))
        
        for GO_id, saliencies in id_saliency_map.items():
            GO_id = GO_id.replace(':', '_') # to stay consistent with other files
            output_file = f"{deepfri_dir}/deepfri_v1_{protein}__{GO_id}.csv"

            # Set residue index to start from 1 instead of 0
            residues = [residue + 1 for residue in range(len(saliencies))]

            df = pd.DataFrame({
                'residue_idx': residues,
                'saliency': saliencies,
            }, index = residues)

            df.to_csv(output_file, index=False)

create_seperate_files(saliency_df, deepfri_dir)
