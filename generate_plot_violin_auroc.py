import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def generate_auroc_violin_plot(input_tsv, output_png, title="AUROC Distributions by Metric"):    
    # 1. Load Data
    if not os.path.exists(input_tsv):
        print(f"Error: File {input_tsv} not found.")
        return
        
    df = pd.read_csv(input_tsv, sep='\t')
    
    # Identify AUROC columns
    metric_cols = [c for c in df.columns if c.endswith('_AUROC')]
    
    # Drop rows with NaN to ensure fair comparison
    df_metrics = df[metric_cols].dropna()
    
    # 2. Sort metrics by their Median AUROC (Best to Worst)
    medians = df_metrics.median().sort_values(ascending=False)
    sorted_metrics = medians.index.tolist()
    
    # 3. Melt the dataframe for Seaborn
    df_melted = df_metrics.melt(var_name='Metric', value_name='AUROC')
    
    # Clean up the metric names for the plot labels (remove '_AUROC' suffix)
    df_melted['Metric'] = df_melted['Metric'].str.replace('_AUROC', '')
    sorted_metrics_clean = [m.replace('_AUROC', '') for m in sorted_metrics]
    
    # Determine Modality for color coding
    def get_modality(metric_name):
        if metric_name.startswith('Seq_'):
            return 'Sequence Model'
        elif metric_name.startswith('Struct_'):
            return 'Structure Model'
        return 'Other'
        
    df_melted['Modality'] = df_melted['Metric'].apply(get_modality)
    
    # 4. Plotting
    # Height is dynamic based on having 20 metrics so it doesn't squish
    plt.figure(figsize=(12, min(12, len(metric_cols) * 0.6)))
    sns.set_theme(style="whitegrid")
    
    ax = sns.violinplot(
        data=df_melted, 
        y='Metric', 
        x='AUROC', 
        hue='Modality',
        dodge=False, # Keeps violins centered on their y-ticks
        order=sorted_metrics_clean,
        palette={'Sequence Model': '#1f77b4', 'Structure Model': '#ff7f0e'}, # Blue & Orange
        inner='quartile', # Draws dashed lines inside the violin for 25%, 50%, 75% quartiles
        cut=0, # Prevents the violin tails from extending past the actual min/max data points
        linewidth=1.2,
        alpha=0.8
    )
    
    # Add a vertical dashed line for random chance
    plt.axvline(x=0.5, color='red', linestyle='--', linewidth=2, alpha=0.8, label='Random Chance (0.5)')
    
    # Formatting
    plt.title(title, fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('AUROC Score', fontsize=14, fontweight='bold')
    plt.ylabel('') # Leave y-label blank since metric names are self-explanatory
    plt.xlim(0, 1.05) # AUROC is strictly 0 to 1
    
    # Tweak the x-axis to show major increments clearly
    plt.xticks([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0], fontsize=11)
    plt.yticks(fontsize=11)
    
    # Handle Legend
    handles, labels = ax.get_legend_handles_labels()
    # Move legend outside the plot so it doesn't cover data
    plt.legend(handles, labels, title='Modality', loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)
    
    plt.tight_layout()
    
    # Save the output
    os.makedirs(os.path.dirname(output_png) or ".", exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Success! Plot saved to: {output_png}\n")

# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Plot Specific Terms
generate_auroc_violin_plot(
    input_tsv="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_specificGOterm_AUROC.tsv",
    output_png="integrated_gradients/auroc_plots/violin_specific_auroc.png",
    title="AUROC Distributions by Metric (Specific GO Terms)"
)

# 2. Plot Propagated Terms
generate_auroc_violin_plot(
    input_tsv="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_propagatedGOterm_AUROC.tsv",
    output_png="integrated_gradients/auroc_plots/violin_propagated_auroc.png",
    title="AUROC Distributions by Metric (Propagated GO Terms)"
)