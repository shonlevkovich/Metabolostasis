# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:25:39 2024

@author: Shon Levkovich
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind, ttest_rel, mannwhitneyu


metabolite_data = pd.read_csv("April combined\onlytestmetabolites_median_norm.csv")  # Replace 'your_file.csv' with the actual file path

# Function to determine the appropriate asterisk annotation for p-values
def p_value_annotation(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"  # not significant

# Melting the data for analysis
metabolite_data_melted = metabolite_data.melt(id_vars=['Metabolite', 'Condition'], 
                                              #value_vars=['1', '2', '3', '4', '5'],
                                              var_name='Repeat', value_name='Measurement')

# Calculating means and standard deviations
metabolite_stats = metabolite_data_melted.groupby(['Metabolite', 'Condition']).agg(['mean', 'std']).reset_index()
metabolite_stats.columns = ['Metabolite', 'Condition', '','', 'Mean', 'Std']

# # Calculate fold changes
# fold_changes = {}
# unique_metabolites = metabolite_data['Metabolite'].unique()
# for metabolite in unique_metabolites:
#     conditions = metabolite_stats[metabolite_stats['Metabolite'] == metabolite]
#     if len(conditions) == 2:
#         mean_control = conditions[conditions['Condition'] == conditions['Condition'].iloc[0]]['Mean'].values[0]
#         mean_test = conditions[conditions['Condition'] == conditions['Condition'].iloc[1]]['Mean'].values[0]
#         fold_change = mean_test / mean_control
#         fold_changes[metabolite] = fold_change
        
# Calculating p-values for the statistical tests
p_values = {}
unique_metabolites = metabolite_data['Metabolite'].unique()
for metabolite in unique_metabolites:
    conditions = metabolite_data[metabolite_data['Metabolite'] == metabolite]['Condition'].unique()
    if len(conditions) == 2:
        data1 = metabolite_data_melted[(metabolite_data_melted['Metabolite'] == metabolite) & 
                                       (metabolite_data_melted['Condition'] == conditions[0])]['Measurement']
        data2 = metabolite_data_melted[(metabolite_data_melted['Metabolite'] == metabolite) & 
                                       (metabolite_data_melted['Condition'] == conditions[1])]['Measurement']
        #t_stat, p_val = ttest_ind(data1, data2, nan_policy='omit')   #, add equal_var=False for Welch t-test
        u_statistic, p_val = mannwhitneyu(data1,data2, nan_policy='omit', alternative='two-sided')
        p_values[metabolite] = p_val

# Plotting
plt.style.use('default')  # Default style for a white background
# Original:
# fig, axes = plt.subplots(1, len(unique_metabolites), figsize=(3 * len(unique_metabolites), 5))

# Modified version for 2 rows:
rows = 2  # Number of rows desired
cols = (len(unique_metabolites) + 1) // rows  # Calculate columns needed, adding 1 to handle odd number of plots properly

fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 5 * rows))  # Adjust figsize accordingly
axes = axes.flatten()  # Flatten the axes array for easy iteration

# Continue with the for loop as before but use the flattened axes
for i, metabolite in enumerate(unique_metabolites):
    ax = axes[i]  # Access the current axis
    subset = metabolite_stats[metabolite_stats['Metabolite'] == metabolite]
    # Assume first condition is control and second is test for fold change calculation
    if len(subset) == 2:
        mean_control = subset[subset['Condition'] == subset['Condition'].iloc[0]]['Mean'].values[0]
        mean_test = subset[subset['Condition'] == subset['Condition'].iloc[1]]['Mean'].values[0]
        fold_change = mean_test / mean_control
        fold_change_text = f'FC={fold_change:.2f}'
    else:
        fold_change_text = 'FC=N/A'  # If not exactly two conditions, display not available
    bar_width = 0.4
    bar_positions = np.arange(len(subset['Condition']))
    for j, condition in enumerate(subset['Condition'].unique()):
        condition_data = subset[subset['Condition'] == condition]
        ax.bar(bar_positions[j], condition_data['Mean'], yerr=condition_data['Std'], width=bar_width, capsize=5, align='center')
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(subset['Condition'], fontsize=14)
    ax.set_title(metabolite, fontsize=14)
    ax.set_ylabel('Normalized concentration (mM/OD600)' if i % cols == 0 else '', fontsize=14)  # Only add ylabel to the first column
    if metabolite == 'Glycine':
        ax.set_ylim(bottom=0)
    p_val = p_values.get(metabolite, 1)
    ax.tick_params(axis='y', labelsize=13)  # Increase y-axis label font size
    ax.text(bar_positions.mean(), max(subset['Mean'] + subset['Std']), p_value_annotation(p_val), 
            ha='center', va='bottom', color='black')
    ax.text(0.5, 0.9, fold_change_text, transform=ax.transAxes, ha='center', color='black', fontsize=10)

# Hide any unused axes if the total number of subplots is less than rows*cols
for i in range(len(unique_metabolites), rows * cols):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig('metabolomics1_april_median_norm.svg', format='svg')
plt.show()


