import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data_file_path = '/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/viral_load_quant/PRRSV2/viral_load_tpm_results_NC134.tsv'
data = pd.read_csv(data_file_path, sep='\t')

# Load the GTF file and create a mapping from gene_id to gene name
gtf_file_path = '/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/ref_genomes/PRRSV2/NC134/NC134_annot.gtf'
gtf_data = pd.read_csv(gtf_file_path, sep='\t', header=None, comment='#',
                       usecols=[0, 2, 8], names=['seqname', 'feature', 'attributes'])
# Extract gene_id and gene name
gtf_data = gtf_data[gtf_data['feature'] == 'gene']
gtf_data['gene_id'] = gtf_data['attributes'].str.extract('gene_id "([^"]+)"')
gtf_data['gene_name'] = gtf_data['attributes'].str.extract('gene "([^"]+)"')
print(gtf_data)
gene_id_to_name = pd.Series(gtf_data['gene_name'].values, index=gtf_data['gene_id']).to_dict()

# Replace gene IDs with gene names in the original data
data['Gene'] = data['Gene'].map(gene_id_to_name)

# Pivot the data for the heatmap
heatmap_data = data.pivot_table(index='Gene', columns='Sample', values='TPM', aggfunc='mean')

# Adjust the data by adding 1 to avoid log(0) which is undefined
adjusted_heatmap_data = np.log1p(heatmap_data)

# Define the custom order for gene names based on gene IDs
gene_id_order = ["ON844087.7", "ON844087.9", "ON844087.11", "ON844087.13", "ON844087.15", "ON844087.17", "ON844087.19", "ON844087.21", "ON844087.23"]
gene_name_order = [gene_id_to_name[gene_id] for gene_id in gene_id_order if gene_id in gene_id_to_name]

# Reorder the heatmap data according to the specified gene order
ordered_heatmap_data = adjusted_heatmap_data.reindex(gene_name_order)

plt.rcParams.update({'font.size': 14})

# Create the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(ordered_heatmap_data, annot=False, cmap='viridis', linewidths=.5, vmin=0, vmax=10, yticklabels=1)
plt.title('Log-Transformed TPM Values by Gene Across PRRSV2 NC134 Samples')
plt.ylabel('Gene')
plt.xlabel('Sample')
plt.yticks(rotation=0)
plt.tight_layout()

# Save the plot as an SVG file
plt.savefig('/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/viral_load_quant/PRRSV2/viral_load_heatmap_NC134_by_gene_name.svg', format='svg')

print("Done.")
