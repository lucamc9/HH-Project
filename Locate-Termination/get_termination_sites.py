'''
Script that given bam file(s) from the command line creates a txt file with
gene names and termination location sites.
s1442231
'''

from joblib import Parallel, delayed
import multiprocessing
from utils import *
import time

# Get bed(s) and bams
bed, utr, name, bams = get_cml_args(sys.argv)

# Safety check
if utr is None:
    print('3\' UTR BED file has not been provided')
    sys.exit()
elif name is None:
    print('No name has been provided for the output')
    sys.exit()
elif bed is None:
    print('BED file has not been provided')
    sys.exit()
elif not bams:
    print('No BAM files have been provided')
    sys.exit()

# Parallelise operations
start_time = time.time()
num_cores = multiprocessing.cpu_count()
term_list = Parallel(n_jobs=num_cores)(delayed(get_termination_df)(bed, utr, bam) for bam in bams)
print('Finish time: {}'.format((time.time()-start_time)/60))

# Create tab delimited file
create_txt_file(name, term_list)
