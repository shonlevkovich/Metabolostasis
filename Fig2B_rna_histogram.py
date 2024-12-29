import pandas as pd
import matplotlib.pyplot as plt

# Load the data file
file_path = 'RNA histogram.xlsx'  # Replace with the correct path to your file
sheet_name = 'Proteomics disribution'
data = pd.read_excel(file_path, sheet_name=sheet_name)

# Extract data for plotting
amino_acids = data['Unnamed: 0']
values_p = data['p<0.05 all']
values_up_down = data['Up/Down >1.5-fold']
values_upregulated = data['Upregulated >1.5-fold']
values_downregulated = data['Downregulated >1.5-fold']

# Custom colors (based on reference image)
colors = {
    'p<0.05 all': '#d9d9d9',  # Light gray
    'Up/Down >1.5-fold': '#6a51a3',  # Purple
    'Upregulated >1.5-fold': '#d7191c',  # Red
    'Downregulated >1.5-fold': '#2c7fb8'  # Blue
}

# Plot parameters
bar_width = 0.2
index = range(len(amino_acids))

# Create the plot
plt.figure(figsize=(10, 5))
plt.bar(index, values_p, bar_width, label='p<0.05 all', color=colors['p<0.05 all'], edgecolor='black')
plt.bar([i + bar_width for i in index], values_up_down, bar_width, label='Up/Down >1.5-fold', color=colors['Up/Down >1.5-fold'], edgecolor='black')
plt.bar([i + 2 * bar_width for i in index], values_upregulated, bar_width, label='Upregulated >1.5-fold', color=colors['Upregulated >1.5-fold'], edgecolor='black')
plt.bar([i + 3 * bar_width for i in index], values_downregulated, bar_width, label='Downregulated >1.5-fold', color=colors['Downregulated >1.5-fold'], edgecolor='black')

# Add labels and legend
plt.xlabel('Amino Acids', fontsize=14)
plt.ylabel('Number of transcripts', fontsize=14)
plt.xticks([i + 1.5 * bar_width for i in index], amino_acids, fontsize=14)
plt.legend(fontsize=14, loc='upper right')

# Remove gridlines and adjust spines
plt.gca().grid(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['bottom'].set_linewidth(1.5)

# Adjust layout and display/save the plot
plt.tight_layout()
plt.savefig('RNA_histogram_publication_grade.svg', dpi=300, bbox_inches='tight')
plt.show()

