# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:25:39 2024

@author: Shon Levkovich
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind, ttest_rel

metabolite_data = pd.read_csv("gptmet.csv")  # Replace 'your_file.csv' with the actual file path

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
                                              value_vars=['1', '2', '3', '4', '5'],
                                              var_name='Repeat', value_name='Measurement')

# Calculating means and standard deviations
metabolite_stats = metabolite_data_melted.groupby(['Metabolite', 'Condition']).agg(['mean', 'std']).reset_index()
metabolite_stats.columns = ['Metabolite', 'Condition', '','', 'Mean', 'Std']

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
        t_stat, p_val = ttest_rel(data1, data2, nan_policy='omit')   ##ttest_rel for paired-samples
        p_values[metabolite] = p_val

# Plotting
plt.style.use('default')  # Default style for a white background
fig, axes = plt.subplots(1, len(unique_metabolites), figsize=(3 * len(unique_metabolites), 5))  # Adjusting subplot size

for i, (metabolite, ax) in enumerate(zip(unique_metabolites, axes), 1):
    subset = metabolite_stats[metabolite_stats['Metabolite'] == metabolite]
    bar_width = 0.3
    bar_positions = np.arange(len(subset['Condition']))
    for j, condition in enumerate(subset['Condition'].unique()):
        condition_data = subset[subset['Condition'] == condition]
        ax.bar(bar_positions[j], condition_data['Mean'], yerr=condition_data['Std'], width=bar_width, capsize=5, align='center')
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(subset['Condition'])
    ax.set_title(metabolite)
    ax.set_ylabel('Intracellular concentration [mM]' if i == 1 else '')
    if metabolite == 'Glycine':
        ax.set_ylim(bottom=0)
    p_val = p_values.get(metabolite, 1)
    ax.text(bar_positions.mean(), max(subset['Mean'] + subset['Std']), p_value_annotation(p_val), 
            ha='center', va='bottom', color='black')

plt.tight_layout()
plt.savefig('metabolomics.svg', format='svg')
plt.show()

