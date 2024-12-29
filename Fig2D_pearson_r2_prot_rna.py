import pandas as pd
import matplotlib.pyplot as plt

# Load the Excel file
file_path = 'r2 diagram.xlsx'  # Replace with the path to your file
sheet_name = 'Proteomics disribution'
data = pd.read_excel(file_path, sheet_name=sheet_name)

# Extract data for plotting
amino_acids = data['Unnamed: 0']
r_all_data = data['r (all data)']
r_significant_only = data['r (significant only)']

# Plot parameters
bar_width = 0.35
index = range(len(amino_acids))

# Create the plot with reduced y-axis aspect ratio
plt.figure(figsize=(8, 4))
plt.bar(index, r_all_data, bar_width, label='r (all data)', color='lightgray', edgecolor='black')
plt.bar([i + bar_width for i in index], r_significant_only, bar_width, label='r (significant only)', color='black', edgecolor='black')

# Add labels, title, and legend
plt.xlabel('Amino Acids', fontsize=10)
plt.ylabel('Pearson correlation coefficient\n(protein vs. mRNA expression change)', fontsize=10)
plt.xticks([i + bar_width / 2 for i in index], amino_acids, fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10, loc='upper left', bbox_to_anchor=(0.01, 1.05))

# Remove gridlines and adjust spines
plt.gca().grid(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['bottom'].set_linewidth(1.5)

# Layout and display
plt.tight_layout()
plt.savefig('pearson_correlation_plot_reduced_aspect.svg', dpi=300, bbox_inches='tight')
plt.show()
