'''
Utils file for finding termination sites.
s1442231
'''

from scipy.stats import poisson
import pysam
import pandas as pd
import numpy as np
import sys
from numba import jit
import time

def locate_termination(samfile, chrom, strand, gene_exons):
    '''
    Returns the window termination location of the gene in gene_exons
    '''
    # Initialize params
    threshold = 0.01
    s_t = 75
    s_w = 100
    l_g = 1000
    n_w = (l_g - s_w)/s_t + 1

    if strand == '-':
        start = gene_exons['LOC1'].iloc[0]
        end = start + 1000
        N = samfile.count(chrom, start=start, end=end)
        N_w = N * s_w / l_g
        m = poisson.ppf(threshold, N_w)

    else:
        end = gene_exons['LOC2'].iloc[-1]
        start = end - 1000
        N = samfile.count(chrom, start=start, end=end)
        N_w = N * s_w / l_g
        m = poisson.ppf(threshold, N_w)

    if m != 0:
        return slide_window(samfile, chrom, start, end, strand, n_w, m, s_w, s_t)
    else:
        return 1

def slide_window(samfile, chrom, start, end, strand, n_w, m, s_w, s_t):
    '''
    Implements the sliding window algorithm, returning l = mid(w_i-1) + mid(w_i) / 2
    '''

    # Init params
    l = 1
    start_wi = 1

    # Move sliding window according to strand
    if strand == '-':
        loc_wi = start
        prev_mid = l + (s_w/2)

        for i in range(int(n_w)):
            count = samfile.count(chrom, start=loc_wi, end=loc_wi+s_w)
            if count <= m:
                prev_mid = start_wi + (s_w/2)
                start_wi += s_t
            else:
                new_mid = start_wi + (s_w/2)
                l = int(np.ceil((prev_mid + new_mid) / 2))
                return l
            loc_wi += s_t
    else:
        loc_wi = end
        prev_mid = l + (s_w/2)

        for i in range(int(n_w)):
            count = samfile.count(chrom, start=loc_wi-s_w, end=loc_wi)
            if count <= m:
                prev_mid = start_wi + (s_w/2)
                start_wi += s_t
            else:
                new_mid = start_wi + (s_w/2)
                l = int(np.ceil((prev_mid + new_mid) / 2))
                return l
            loc_wi -= s_t

    return l

def get_gene_info(bed, utr, chrom, next_loc):
    '''
    Returns the gene_name and gene_exons for a gene in the region where next_loc is
    '''

    try:
        gene = bed[(bed['CHROM'] == chrom) & (bed['LOC1'] <= next_loc) & (bed['LOC2'] >= next_loc)]
        gene_name = gene['NAME'].iloc[0]
        return gene_name, utr[utr['EXON'].str.contains(gene_name)]
    except IndexError:
        return None, None

def get_gene_end(bed, loc, chrom):
    '''
    Returns the location of the end of the gene given location loc and chromosome chrom
    '''

    gene = bed[(bed['CHROM'] == chrom) & (bed['LOC1'] <= loc) & (bed['LOC2'] >= loc)]
    return gene['LOC2'].iloc[0]

@jit
def get_termination_df(bed_name, utr_name, bam_name):
    '''
    Returns a DataFrame with all columns of argument bed plus a column of termination locations
    '''

    # Load bed files
    bed = pd.read_csv(bed_name, delimiter='\t',
                    names=['CHROM', 'LOC1', 'LOC2', 'NAME', 'STRAND'],
                    converters={'CHROM': str.strip, 'NAME': str.strip, 'STRAND': str.strip})
    utr = pd.read_csv(utr_name, delimiter='\t',
                    names=['CHROM', 'LOC1', 'LOC2', 'EXON', 'STRAND'])

    # Initialize params
    samfile = pysam.AlignmentFile(bam_name, "rb")
    term_bed = bed
    term_bed['TERM'] = 0
    chromosomes = ['1','2','3','4','5','6','7','8','9',
                   '10','11','12','13','14','15','16',
                   '17', '18', '19', 'X', 'Y']

    # Iterate over chromosomes and reads
    start_time = time.time()
    for chrom in chromosomes:
        border = 0
        for read in samfile.fetch(chrom):
            start_position = read.get_reference_positions()[0]
            if start_position <= border: # Skip read if in region already accounted for
                continue
            gene_name, gene_exons = get_gene_info(bed, utr, chrom, start_position)
            if gene_exons is not None and not gene_exons.empty:
                strand = gene_exons['STRAND'].iloc[0]
                term_location = locate_termination(samfile, chrom, strand, gene_exons)
                term_bed.loc[term_bed['NAME'] == gene_name, 'TERM'] = term_location
                border = get_gene_end(bed, start_position, chrom)
            else:
                continue

    print('{} ...done'.format(bam_name))
    print('Time taken: {} mins'.format((time.time() - start_time)/60))
    return term_bed

def get_cml_args(cml_args):
    '''
    Command line parser
    '''

    bed = None
    utr = None
    name = None
    bams = []
    for i in range(1, len(cml_args)):
        indicator = cml_args[i-1]
        if indicator == '--bam':
            bams = cml_args[i:len(cml_args)]
            break
        elif indicator == '--utr':
            utr = cml_args[i]
        elif indicator == '--bed':
            bed = cml_args[i]
        elif indicator == '--name':
            name = cml_args[i]

    return bed, utr, name, bams

def create_txt_file(name, term_list):
    '''
    Creates txt files with the gene names as the first column and termination location
    of bams following
    '''

    df = term_list[0][['NAME', 'TERM']]
    for i in range(1, len(term_list)):
        df[str(i)] = term_list[i][['TERM']]

    df.to_csv(name +'.txt', header=None, index=None, sep='\t', mode='a')
