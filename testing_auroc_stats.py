import pandas as pd
import scipy.stats as stats
import scikit_posthocs as sp
import matplotlib.pyplot as plt
import os

def run_statistical_analysis(auroc_file, output_plot_path):
    
    # 1. Load Data
    if not os.path.exists(auroc_file):
        print(f"Error: File {auroc_file} not found.")
        return
        
    df = pd.read_csv(auroc_file, sep='\t')
    
    # Isolate only the AUROC columns
    metric_cols = [col for col in df.columns if col.endswith('_AUROC')]
    
    # The Friedman test requires "complete blocks" (no missing data)
    df_metrics = df[metric_cols].dropna()
    
    print(f"Analyzing {len(df_metrics)} valid protein targets across {len(metric_cols)} metrics.")
    
    # 2. The Friedman Test
    # We unpack the dataframe columns into separate arguments for the SciPy function
    stat, p_value = stats.friedmanchisquare(*[df_metrics[col] for col in metric_cols])
    print(f"Friedman Test Results: statistic = {stat:.3f}, p-value = {p_value:.3e}")
    
    if p_value < 0.05:
        print("Result: Significant differences exist between the metrics. Proceeding to post-hoc analysis...")
        
        # 3. Calculate Average Ranks
        # We rank across axis=1 (across the metrics for each protein). 
        # ascending=False means the highest AUROC gets Rank 1.
        ranks = df_metrics.rank(axis=1, ascending=False)
        avg_ranks = ranks.mean()
        
        # 4. The Nemenyi Post-Hoc Test
        # This calculates the p-value matrix comparing every metric to every other metric
        nemenyi_p_values = sp.posthoc_nemenyi_friedman(df_metrics.values)
        nemenyi_p_values.columns = metric_cols
        nemenyi_p_values.index = metric_cols
        
        # 5. Generate the Critical Difference (CD) Diagram
        plt.figure(figsize=(10, 6))
        
        # Clean up column names for the plot labels (remove '_AUROC' for readability)
        clean_labels = [col.replace('_AUROC', '') for col in metric_cols]
        avg_ranks.index = clean_labels
        nemenyi_p_values.columns = clean_labels
        nemenyi_p_values.index = clean_labels
        
        # Draw the diagram
        sp.critical_difference_diagram(
            ranks=avg_ranks, 
            sig_matrix=nemenyi_p_values,
            label_fmt_left='{label} [{rank:.2f}]', 
            label_fmt_right='[{rank:.2f}] {label}',
            text_h_margin=0.2  # Gives labels some breathing room
        )
        
        plt.title(f"Critical Difference Diagram\n(Nemenyi test, $\\alpha=0.05$)", pad=20, fontweight='bold')
        
        # Save the plot
        os.makedirs(os.path.dirname(output_plot_path) or ".", exist_ok=True)
        plt.savefig(output_plot_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"Success! Critical Difference diagram saved to: {output_plot_path}")
        
    else:
        print("Result: No statistically significant differences found between the metrics. (p >= 0.05)")

# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Analyze Specific Terms
run_statistical_analysis(
    auroc_file="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_specificGOterm_AUROC.tsv",
    output_plot_path="integrated_gradients/auroc_plots/cd_diagram_specific.png"
)

# 2. Analyze All/Propagated Terms
run_statistical_analysis(
    auroc_file="integrated_gradients/biolip_multimetric_overlap_combined_residues_ligand_propagatedGOterm_AUROC.tsv",
    output_plot_path="integrated_gradients/auroc_plots/cd_diagram_propagated.png"
)
