'''
Script that merges a gene expression matrix (GEM) and a gene termination matrix (GTM)
into a merged GTM, either as tupled values or double rows
s1442231
'''

import pandas as pd
import numpy as np
from numba import jit
import sys

# Getting command line args
cml_args = sys.argv[1:]

gtm = None
gem = None
as_tuple = False
filter_genes = False
name = None
X = 0

for i in range(1, len(cml_args)):
    indicator = cml_args[i-1]
    if indicator == '--gtm':
        gtm = cml_args[i]
    elif indicator == '--gem':
        gem = cml_args[i]
    elif indicator == '--tuple':
        as_tuple = True
    elif indicator == '--name':
        name = cml_args[i]
    elif indicator == '--filter':
        filter_genes = True
        X = int(cml_args[i])

# Safety check
if gem is None:
    print('No GEM has been provided')
    sys.exit()
elif name is None:
    print('No name has been provided for the output')
    sys.exit()
elif gtm is None:
    print('No GTM has been provided')
    sys.exit()


# Initialize dataframes
gem = pd.read_csv(gem, delimiter='\t', header=None)
gtm = pd.read_csv(gtm, delimiter='\t', header=None)

# Remove artefact
gem = gem[gem[0] != 'AK040785']
gtm = gtm.drop_duplicates()

if as_tuple:
    # Merge & make sure there's no NaN values
    merged_df = pd.merge(gem, gtm, on=0, how='left').dropna(axis=0, how='any')

    # Separate sorted matrices
    term_df = merged_df.iloc[:, gem.shape[1]:]
    expr_df = merged_df.iloc[:, 1:gem.shape[1]]

    final_df = merged_df[[0]]   # Initialize df with gene_names
    for i in range(len(expr_df.columns)):
        final_df[i+1] = list(zip(expr_df.iloc[:,i], term_df.iloc[:,i].astype(int)))
    name = 'tpl_' + name

else:
    # Sort gtm as gem
    gtm = gtm.set_index(0)
    gem_copy = gem.copy()
    gem_copy = gem_copy.set_index(0)
    gtm = gtm.reindex(gem_copy.index).dropna(axis=0, how='any').astype(int)

    # Make sure gtm values aligns properly with gem
    bool_matrix = gem_copy.apply(lambda x: x > 0)
    gtm = gtm * bool_matrix

    # Filter rare genes
    if filter_genes:
        num_cells = gem.shape[1]
        X = X/100
        filter_ = X*num_cells
        # condition
        bool_expr = gem_copy.apply(lambda x: x > 2)
        # sum of true values
        bool_expr = bool_expr.sum(axis=1)
        # apply filter
        genes_kept = bool_expr[bool_expr > filter_].index
        filtered_gem = gem_copy.loc[genes_kept, :]
        filtered_gtm = gtm.loc[genes_kept,:]
        # Change index back to numerical
        gtm = filtered_gtm.reset_index()
        gem = filtered_gem.reset_index()
    else:
        # Change index back to numerical
        gtm = gtm.reset_index()

    # Append_e or _t to all values in first column
    gem[0] = gem[0].astype(str) + '_e'
    gtm[0] = gtm[0].astype(str) + '_t'

    # Set indexes even/odd
    idx_even = list(range(0, gem.shape[0]*2 - 1, 2))
    idx_odd = list(range(1, gem.shape[0]*2, 2))
    gem = gem.set_index([idx_even])
    gtm = gtm.set_index([idx_odd])

    # Concatenate dfs and sort index
    final_df = pd.concat([gem, gtm]).sort_index()

# Write to txt
final_df.to_csv(name, header=None, index=None, sep='\t', mode='a')
