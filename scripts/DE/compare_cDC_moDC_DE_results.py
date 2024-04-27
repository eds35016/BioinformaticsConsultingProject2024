import pandas as pd

# Load the TSV files
cdc_data = pd.read_csv('/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/DE/STAR_counts/PRRSV2_adjusted_postQC/cDC_de.m.tsv', delimiter='\t')
modc_data = pd.read_csv('/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/DE/STAR_counts/PRRSV2_adjusted_postQC/moDC_de.m.tsv', delimiter='\t')

# Extracting significantly differentially expressed genes for each comparison in both datasets
# cDC data
cdc_134_vs_mock = set(cdc_data[cdc_data['cDC_134.vs.Mock.flag']].index)
cdc_174_vs_mock = set(cdc_data[cdc_data['cDC_174.vs.Mock.flag']].index)
cdc_134vmock_vs_174vmock = set(cdc_data[cdc_data['cDC.134vMock.vs.174vMock.flag']].index)

# moDC data
modc_134_vs_mock = set(modc_data[modc_data['moDC_134.vs.Mock.flag']].index)
modc_174_vs_mock = set(modc_data[modc_data['moDC_174.vs.Mock.flag']].index)
modc_134vmock_vs_174vmock = set(modc_data[modc_data['moDC.134vMock.vs.174vMock.flag']].index)

# Calculating overlaps
overlap_134_vs_mock = len(cdc_134_vs_mock & modc_134_vs_mock)
overlap_174_vs_mock = len(cdc_174_vs_mock & modc_174_vs_mock)
overlap_134vmock_vs_174vmock = len(cdc_134vmock_vs_174vmock & modc_134vmock_vs_174vmock)

# Output the results
print(f'Overlap in 134 vs Mock: {overlap_134_vs_mock}')
print(f'Overlap in 174 vs Mock: {overlap_174_vs_mock}')
print(f'Overlap in 134vMock vs 174vMock: {overlap_134vmock_vs_174vmock}')
