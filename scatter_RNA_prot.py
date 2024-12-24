# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:35:51 2024

@author: Shon Levkovich
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np, math
from sklearn.linear_model import LinearRegression
from adjustText import adjust_text
from matplotlib_venn import venn2
from math import log2
import os


# Load the data
file_path = 'C:/Users/Shon Levkovich/OneDrive - Tel-Aviv University/TAU/RNA-Proteomics/3 - Joint analysis proteomic-RNA (scatter and regression)/1 - 2D plots/Input/Met.xlsx'
data = pd.read_excel(file_path)

# Extract the base name (file name with extension)
base_name = os.path.basename(file_path)

# Remove the file extension to get just the name
metabolite = os.path.splitext(base_name)[0]


# Clean and prepare data
data.columns = ['Gene name Proteomics', 'Proteomics p (-log10)', 'Proteomics FC (log2)', 
                'Gene name Transcriptomics', 'Transcriptomics p (-log10)', 'Transcriptomics FC (log2)']
data['Proteomics p (-log10)'] = pd.to_numeric(data['Proteomics p (-log10)'], errors='coerce')
data['Proteomics FC (log2)'] = pd.to_numeric(data['Proteomics FC (log2)'], errors='coerce')
data['Transcriptomics p (-log10)'] = pd.to_numeric(data['Transcriptomics p (-log10)'], errors='coerce')
data['Transcriptomics FC (log2)'] = pd.to_numeric(data['Transcriptomics FC (log2)'], errors='coerce')

# Merge on gene name
merged_data = pd.merge(data[['Gene name Proteomics', 'Proteomics p (-log10)', 'Proteomics FC (log2)']], 
                       data[['Gene name Transcriptomics', 'Transcriptomics p (-log10)', 'Transcriptomics FC (log2)']], 
                       left_on='Gene name Proteomics', right_on='Gene name Transcriptomics', 
                       how='inner')

# Calculate p-values from -log10 p-values
merged_data['Proteomics_p_value'] = 10**(-merged_data['Proteomics p (-log10)'])
merged_data['Transcriptomics_p_value'] = 10**(-merged_data['Transcriptomics p (-log10)'])

# Determine significance
threshold = -math.log10(0.05)  # Corresponds to p-value ~ 0.05
merged_data['Significance'] = (
    (merged_data['Proteomics_p_value'] <= 0.05).astype(int) + 
    (merged_data['Transcriptomics_p_value'] <= 0.05).astype(int) * 2
)

# Color mapping
color_map = {0: 'Not significant', 1: 'Significant only in proteomics', 2: 'Significant only in transcriptomics', 3: 'Significant in both'}
merged_data['Color'] = merged_data['Significance'].map(color_map)

# Drop any rows with NaN values in columns used for regression
cleaned_data = merged_data.dropna(subset=['Proteomics FC (log2)', 'Transcriptomics FC (log2)'])

# Regression analysis for all data
x = cleaned_data['Transcriptomics FC (log2)'].values.reshape(-1, 1)
y = cleaned_data['Proteomics FC (log2)'].values
model_all = LinearRegression()
model_all.fit(x, y)
y_pred_all = model_all.predict(x)

# Regression analysis for only the purple data significant in both
sig_both = cleaned_data[cleaned_data['Color'] == 'Significant in both']
x_sig_both = sig_both['Transcriptomics FC (log2)'].values.reshape(-1, 1)
y_sig_both = sig_both['Proteomics FC (log2)'].values
model_sig_both = LinearRegression()
model_sig_both.fit(x_sig_both, y_sig_both)
y_pred_sig_both = model_sig_both.predict(x_sig_both)

# Add a column for the Manhattan distance from the origin to identify the top genes
merged_data['Manhattan Distance'] = merged_data['Proteomics FC (log2)'].abs() + merged_data['Transcriptomics FC (log2)'].abs()

# Identify the top 15 genes based on Manhattan distance
top_genes = merged_data[merged_data['Color'] == 'Significant in both'].nlargest(15, 'Manhattan Distance')

# Plotting
plt.figure(figsize=(10, 8))
colors = {'Not significant': '#cecece', 'Significant only in proteomics': '#ae282c',
          'Significant only in transcriptomics': '#2066a8', 'Significant in both': '#8d3781'}

for key, group in cleaned_data.groupby('Color'):
    plt.scatter(group['Transcriptomics FC (log2)'], group['Proteomics FC (log2)'], label=key, color=colors[key], alpha=0.4)


#draw regression lines
plt.plot(x, y_pred_all, color='black', linestyle='dotted', label='Regression All Data')
plt.plot(x_sig_both, y_pred_sig_both, color='#8d3781', linestyle='dotted', label='Regression Significant Both')

# Calculate R² values and display them on the plot
r2_all = model_all.score(x, y)
r2_sig_both = model_sig_both.score(x_sig_both, y_sig_both)

# Calculate Pearson correlation coefficients for all data
correlation_all = cleaned_data['Proteomics FC (log2)'].corr(cleaned_data['Transcriptomics FC (log2)'])

# Calculate Pearson correlation coefficients for "significant in all three" subset
subset_both = cleaned_data[cleaned_data['Color'] == 'Significant in both']
correlation_sigboth = subset_both['Proteomics FC (log2)'].corr(subset_both['Transcriptomics FC (log2)'])


# Add text annotations for R² values
plt.text(0.05, 0.95, f'R² (all data) = {r2_all:.2f}', transform=plt.gca().transAxes, ha='left', color='black', fontweight='bold')
plt.text(0.05, 0.92, f'Pearson r (all data) = {correlation_all:.2f}', transform=plt.gca().transAxes, ha='left', color='black', fontweight='bold')
plt.text(0.05, 0.89, f'R² (significant both)= {r2_sig_both:.2f}', transform=plt.gca().transAxes, ha='left', color='#8d3781', fontweight='bold')
plt.text(0.05, 0.86, f'Pearson r (significant both) = {correlation_sigboth:.2f}', transform=plt.gca().transAxes, ha='left', color='#8d3781', fontweight='bold')

##Optional - to calculate the Pearson correlation coefficient (r), which for linear regression equals the sqaure root of R^2. 
##Alternatively, you can just take the SQRT of R^2. 
# Calculate Pearson correlation coefficient
#pearson_corr_all = np.corrcoef(cleaned_data['Proteomics FC (log2)'], cleaned_data['Transcriptomics FC (log2)'])[0, 1]
#pearson_corr_sig_both = np.corrcoef(sig_both['Proteomics FC (log2)'], sig_both['Transcriptomics FC (log2)'])[0, 1]

# Add text annotations for Pearson correlation
#plt.text(0.05, 0.85, f'All Data Pearson r = {pearson_corr_all:.2f}', transform=plt.gca().transAxes, ha='left', color='black', fontweight='bold')
#plt.text(0.05, 0.80, f'Significant Both Pearson r = {pearson_corr_sig_both:.2f}', transform=plt.gca().transAxes, ha='left', color='#8d3781', fontweight='bold')

# Plot enhancements
plt.title(f'{metabolite} Proteomics and Transcriptomics Correlation')
plt.xlabel('Log2 Fold Change - Transcriptomics')
plt.ylabel('Log2 Fold Change - Proteomics')
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.legend(loc='center right', bbox_to_anchor=(1.4, 0.5))

texts = []  # List to hold all the text objects for adjust_text
for index, row in top_genes.iterrows():
    plt.scatter(row['Proteomics FC (log2)'], row['Transcriptomics FC (log2)'], color=colors[row['Color']], edgecolors='black', linewidths=0.2, s=40)
    # Create a text object for each gene and append to list
    text = plt.text(row['Proteomics FC (log2)'], row['Transcriptomics FC (log2)'], ' ' + row['Gene name Proteomics'], fontsize=9, ha='right', va='center')
    texts.append(text)

# Use adjust_text to optimize text placement
adjust_text(texts, 
            x=top_genes['Transcriptomics FC (log2)'].values, 
            y=top_genes['Proteomics FC (log2)'].values, 
            force_text=(0.1,0.2),  # Correct usage: a single float value
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))  # Ensure arrows are visible

# Save plot and data
plt.savefig(f'{metabolite}_scatter_plot.svg', format='svg')
cleaned_data.to_csv(f'{metabolite}_scatter_data.csv', index=False)

print('Plot and data saved successfully.')

############################################3
##############Venn plotter
#############################################
# Filter for genes that meet the criteria for proteomics and transcriptomics
fc_cutoff_linear = 1.5
fc_cutoff_log2 = log2(fc_cutoff_linear)

proteomics_criteria = (cleaned_data['Proteomics_p_value'] <= 0.05) & (abs(cleaned_data['Proteomics FC (log2)']) >= fc_cutoff_log2)
transcriptomics_criteria = (cleaned_data['Transcriptomics_p_value'] <= 0.05) & (abs(cleaned_data['Transcriptomics FC (log2)']) >= fc_cutoff_log2)

proteomics_genes = set(cleaned_data[proteomics_criteria]['Gene name Proteomics'])
transcriptomics_genes = set(cleaned_data[transcriptomics_criteria]['Gene name Transcriptomics'])

# Venn Diagram
plt.figure(figsize=(8, 8))
venn_diagram=venn2([proteomics_genes, transcriptomics_genes], ('Proteomics', 'Transcriptomics'),
                   set_colors=('#ae282c', '#2066a8'), alpha=0.5)

# Set label sizes
venn_diagram.get_label_by_id('10').set_fontsize(16)  # Proteomics only
venn_diagram.get_label_by_id('01').set_fontsize(16)  # Transcriptomics only
venn_diagram.get_label_by_id('11').set_fontsize(16)  # Intersection

# Set the labels for the sets and increase their size
for text in venn_diagram.set_labels:
    text.set_fontsize(18)

# Optionally, increase annotation size if needed
for text in venn_diagram.subset_labels:
    if text:  # Check if the label exists (i.e., not None)
        text.set_fontsize(16)
        
plt.title(f'{metabolite} Venn Diagram of Significant Genes (log\u2082|FC|>log\u2082(1.5))', fontsize=18)
plt.savefig(f'{metabolite}_venn.svg', format='svg')


# Collect data for CSV export
overlap_genes = proteomics_genes & transcriptomics_genes
only_proteomics_genes = proteomics_genes - transcriptomics_genes
only_transcriptomics_genes = transcriptomics_genes - proteomics_genes

# Create CSV files with additional FC and p-value data
def extract_gene_data(gene_set, data):
    # Filter data for specific gene set
    gene_data = data[data['Gene name Proteomics'].isin(gene_set)]
    return gene_data[['Gene name Proteomics', 'Proteomics FC (log2)', 'Proteomics_p_value', 'Proteomics p (-log10)', 'Transcriptomics FC (log2)', 'Transcriptomics_p_value', 'Transcriptomics p (-log10)']]

overlap_data = extract_gene_data(overlap_genes, cleaned_data)
proteomics_data = extract_gene_data(only_proteomics_genes, cleaned_data)
transcriptomics_data = extract_gene_data(only_transcriptomics_genes, cleaned_data)

# Save to CSV
overlap_data.to_csv(f'{metabolite}_overlap_genes.csv', index=False)
proteomics_data.to_csv(f'{metabolite}_only_proteomics_genes.csv', index=False)
transcriptomics_data.to_csv(f'{metabolite}_only_transcriptomics_genes.csv', index=False)

###########################################
