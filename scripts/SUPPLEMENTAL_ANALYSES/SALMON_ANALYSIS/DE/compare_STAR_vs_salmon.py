import pandas as pd
import os

# Set paths
STAR_path = '/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/DE/STAR_counts/PRRSV2_adjusted'
salmon_path = '/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/DE/salmon_counts/PRRSV2_adjusted'

# Load the TSV files
df1 = pd.read_csv(os.path.join(STAR_path, 'cDC_sig_counts.tsv'), sep='\t', index_col=0)
df2 = pd.read_csv(os.path.join(salmon_path, 'cDC_sig_counts.tsv'), sep='\t', index_col=0)
df3 = pd.read_csv(os.path.join(STAR_path, 'moDC_sig_counts.tsv'), sep='\t', index_col=0)
df4 = pd.read_csv(os.path.join(salmon_path, 'moDC_sig_counts.tsv'), sep='\t', index_col=0)

def compare_methods(df_a, df_b):
    """Compare differential expression methods by calculating Pearson correlation for common genes."""
    # Ensure the samples are in the same order for comparison
    df_a = df_a.reindex(sorted(df_a.columns), axis=1)
    df_b = df_b.reindex(sorted(df_b.columns), axis=1)
    
    # Check if the number of samples is the same
    if df_a.shape[1] != df_b.shape[1]:
        raise ValueError("The number of samples in each method must be the same for direct comparison.")
    
    common_genes = df_a.index.intersection(df_b.index)
    
    correlations = []
    for col in df_a.columns:
        corr = df_a[col].loc[common_genes].corr(df_b[col].loc[common_genes])
        correlations.append(corr)
        
    # Calculate the mean correlation across all samples
    mean_correlation = sum(correlations) / len(correlations)
    return mean_correlation, common_genes

def calculate_jaccard_index(set_a, set_b):
    """Calculate the Jaccard Index and print set statistics."""
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    jaccard_index = intersection / union
    print(f"Number of significantly DE genes in STAR results: {len(set_a)}")
    print(f"Number of significantly DE genes in salmon results: {len(set_b)}")
    print(f"Intersection size: {intersection}")
    print(f"Union size: {union}")
    return jaccard_index

# Compare the methods for each cell type
correlation_ct1, common_genes_ct1 = compare_methods(df1, df2)
correlation_ct2, common_genes_ct2 = compare_methods(df3, df4)

print("Average logFC correlation between methods for cDC results:", correlation_ct1)
print("Average logFC correlation between methods for moDC results:", correlation_ct2)

# Convert indices to sets for Jaccard Index calculation
genes_set1 = set(df1.index)
genes_set2 = set(df2.index)
genes_set3 = set(df3.index)
genes_set4 = set(df4.index)

# Calculate and print Jaccard Index for each cell type
print("\nJaccard Index calculation for cDC results:")
jaccard_index_ct1 = calculate_jaccard_index(genes_set1, genes_set2)
print("Jaccard Index for cDC:", jaccard_index_ct1)

print("\nJaccard Index calculation for moDC results:")
jaccard_index_ct2 = calculate_jaccard_index(genes_set3, genes_set4)
print("Jaccard Index for moDC:", jaccard_index_ct2)
