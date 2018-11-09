'''
Script that merges gene transcription matrices (GTMs) in order
s1442231
'''

import pandas as pd
import numpy as np
import sys

gtms = sys.argv[1:]

if len(gtms) < 2:
    print('Need more than 1 GTM to merge')
    sys.exit()

merged_gtm = pd.read_csv(gtms[0], delimiter='\t', header=None)
for gtm_name in gtms[1:]:
    gtm = pd.read_csv(gtm_name, delimiter='\t', header=None)
    merged_gtm = pd.concat([merged_gtm, gtm.iloc[:, 1:]], axis=1, ignore_index=True)

merged_gtm.to_csv('gtm_final.txt', header=None, index=None, sep='\t', mode='a')
