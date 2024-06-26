import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
file_path = '/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/viral_load_quant/PRRSV2/viral_load_tpm_results_NC174.tsv'
data = pd.read_csv(file_path, sep='\t')

# Pivot the data for the heatmap
heatmap_data = data.pivot_table(index='Gene', columns='Sample', values='TPM', aggfunc='mean')

# Adjust the data by adding 1 to avoid log(0) which is undefined
adjusted_heatmap_data = np.log1p(heatmap_data)

# Define the custom order for genes
gene_order = ["NC18-9-7.6", "NC18-9-7.8", "NC18-9-7.10", "NC18-9-7.12", "NC18-9-7.14", "NC18-9-7.16", "NC18-9-7.18"]

# Reorder the heatmap data according to the specified gene order
ordered_heatmap_data = adjusted_heatmap_data.reindex(gene_order)

plt.rcParams.update({'font.size': 14})

# Create the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(ordered_heatmap_data, annot=False, cmap='viridis', linewidths=.5, vmin=0, vmax=10)
plt.title('Log-Transformed TPM Values by Gene Across PRRSV2 NC174 Samples')
plt.ylabel('Gene')
plt.xlabel('Sample')
plt.tight_layout()

# Save the plot as an SVG file
plt.savefig('/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/viral_load_quant/PRRSV2/viral_load_heatmap_NC174.svg', format='svg')

print("Done.")
