## Analyze Fold Switching Proteins

#### Inputs
- `raw_data`: predictions and scores for two folds of the same protein from DeepFRI v1.0 and DeepFRI v2.0 (separated by model)
- `fold_switching_metadata.tsv`: protein names, fold names, metadata *(e.g. sequence identity, TM-score)


#### Reproduce
1. Run analysis.py to reproduce `results`
2. Run `plot_single_protein.py` to reproduce CopK in-depth analysis
3. Run `scatter_plots.py` to reproduce plots in `plots`
4. Run `IC_comparison/workflow.sh` to reproduce semantic similarity comparison **(requires Groovy)**


#### Sub-directory structure:
- `utils`: sym link to semantic utilities
- `raw_data`:
    - deepfriv1_data: raw output from DeepFRI v1.0
    - models_propagated_data: raw output from DeepFRI v2.0
- `plots`: plots generated from 'plot_*' scripts
- `results`:
    - results_fold_switching_deltas.tsv: only terms that are above the threshold in either fold1 or fold2, delta is difference in score (score fold1 - score fold2)
    - results_fold_switching_jaccard.tsv: Jaccard score for the set of terms that pass threshold for each model (comparing folds)
    - results_fold_switching_threshold_crossings.tsv: terms that move below or above the threshold for each model depending on the fold
    - results_fold_switching_semantic_sim.tsv: Lin's semantic similarity + Best Match Average
- `IC_comparison`: semantic similarity score calculation ([Resnik or Lin] + Best Match Average (BMA))
