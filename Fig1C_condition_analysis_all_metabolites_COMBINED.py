# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:53:32 2024

@author: Shon Levkovich
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# Load your dataset
file_path = 'Metaboanalyst/Normalization by median and log10, fill missing with median/1 - Normalization results/data_normalized_only_median.csv'  # Update this to your dataset file path
data = pd.read_csv(file_path)
# Pre-process data: Exclude specified metabolites and handle missing values as needed
data = data.drop(columns=['HO-Lys', 'b-Ala'])



# Define your conditions and metabolites
conditions = data['Label'].unique()
conditions = [cond for cond in conditions if 'WT+' in cond]
metabolites = data.columns[2:]


# Analyze data: Calculate means, STDs, and perform t-tests
analysis_results = {}
for condition in conditions:
    condition_data = {}
    for metabolite in metabolites:
        control_values = data[data['Label'] == 'WT'][metabolite].values
        treatment_values = data[data['Label'] == condition][metabolite].values
        if len(control_values) > 1 and len(treatment_values) > 1:
            control_mean, treatment_mean = np.mean(control_values), np.mean(treatment_values)
            control_std, treatment_std = np.std(control_values, ddof=1), np.std(treatment_values, ddof=1)
            #t_stat, p_value = ttest_ind(control_values, treatment_values) #t_test assumes normal distribution
            u_statistic, p_value = mannwhitneyu(control_values, treatment_values, nan_policy='omit', alternative='two-sided')
            condition_data[metabolite] = {'control_mean': control_mean, 'treatment_mean': treatment_mean,
                                          'control_std': control_std, 'treatment_std': treatment_std,
                                          'p_value': p_value}
    analysis_results[condition] = condition_data
        # At this point in the code, all p_values have been calculated and stored in analysis_results
    all_p_values = []
    for condition_data in analysis_results.values():
        for metabolite_data in condition_data.values():
            all_p_values.append(metabolite_data['p_value'])
    
    # Apply the Benjamini-Hochberg correction
    _, pvals_corrected, _, _ = multipletests(all_p_values, alpha=0.05, method='fdr_bh')
    
    # Update the analysis_results with the corrected p-values
    p_value_index = 0
    for condition_data in analysis_results.values():
        for metabolite_data in condition_data.values():
            metabolite_data['p_value_corrected'] = pvals_corrected[p_value_index]
            p_value_index += 1
    

# Function to add p-value annotations
def add_p_value_annotations(ax, bars, p_values):
    for bar, p in zip(bars, p_values):
        height = bar.get_height()
        p_text = '****' if p<0.0001 else '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        ax.text(bar.get_x() + bar.get_width() / 2., height*1.01, p_text, ha='center', va='bottom', fontsize=10)

# Generate figures with enhancements
plt.rcParams.update({'font.size': 12})  # Adjust font size
fig, axs = plt.subplots(len(conditions), 1, figsize=(8, 2*len(conditions)), dpi=300)  # Adjust for high resolution
fig.subplots_adjust(hspace=0.5)
fig.supylabel('log$_{10}$ concentration (mM)')

for ax, condition in zip(axs, conditions):
    metabolites_list = list(analysis_results[condition].keys())
    control_means = [analysis_results[condition][met]['control_mean'] for met in metabolites_list]
    treatment_means = [analysis_results[condition][met]['treatment_mean'] for met in metabolites_list]
    control_stds = [analysis_results[condition][met]['control_std'] for met in metabolites_list]
    treatment_stds = [analysis_results[condition][met]['treatment_std'] for met in metabolites_list]
    p_values = [analysis_results[condition][met]['p_value_corrected'] for met in metabolites_list]

    bar_width = 0.35
    r1 = np.arange(len(control_means))
    r2 = [x + bar_width for x in r1]

    bars1 = ax.bar(r1, control_means, width=bar_width, yerr=control_stds, capsize=5, label='Control', color='skyblue', edgecolor='black')
    bars2 = ax.bar(r2, treatment_means, width=bar_width, yerr=treatment_stds, capsize=5, label=condition.split('+')[1], color='orange', edgecolor='black')

    add_p_value_annotations(ax, bars2, p_values)

    ax.set_yscale('log')

    ax.set_title(condition)
    ax.set_xticks([r + bar_width/2 for r in range(len(control_means))])
    ax.set_xticklabels(metabolites_list, rotation=45, ha='right')
    # Position legend on the right side outside the subplot
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
    ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig(fname='figure.svg', format='svg')
plt.show()
# Compile analysis results into a DataFrame for saving
results_list = []
for condition, metabolites_data in analysis_results.items():
    for metabolite, stats in metabolites_data.items():
        results_list.append({
            'Condition': condition,
            'Metabolite': metabolite,
            'Control Mean': stats['control_mean'],
            'Treatment Mean': stats['treatment_mean'],
            'Control STD': stats['control_std'],
            'Treatment STD': stats['treatment_std'],
            'P-Value': stats['p_value'],
            'P-Value_BH': stats['p_value_corrected']
        })

results_df = pd.DataFrame(results_list)

# Save the DataFrame to a CSV file
results_csv_path = 'metabolomics_analysis_results.csv'  # Define your desired output path
results_df.to_csv(results_csv_path, index=False)

print(f"Analysis results saved to {results_csv_path}")
