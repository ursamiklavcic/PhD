#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
"""A script to calculate TPM values for contigs or genes based on count files

TPM values are defined as in Wagner et al (Theory in Biosciences) 2012.

      rg x rl x 10^6
TPM = --------------
        flg x T

rg: reads mapped to gene g
rl: read length
flg: feature length
T: sum of rgxrl/flg for all genes
"""

import sys
import pandas as pd
import argparse
import logging

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def main(args):
    logging.info("Reading sample info")
    sample_info = pd.read_table(
        args.sample_info,
        header=None,
        index_col=0,
        names=['avg_read_len']
    )

    logging.info("Reading gene lengths")
    gene_lengths = pd.read_table(
        args.gene_lengths,
        header=None,
        index_col=0,
        names=['gene_length']
    )

    df = pd.DataFrame()

    for fn, sample_name in zip(args.coverage_files, args.sample_names):
        logging.info("Calculating TPM for " + sample_name)

        # Read counts
        rg = pd.read_table(fn, index_col=0, header=None, names=['count'])

        # Keep only genes with length info
        common = rg.index.intersection(gene_lengths.index)
        rg = rg.loc[common]
        gl = gene_lengths.loc[common]

        # Read length for sample
        rl = sample_info.loc[sample_name, 'avg_read_len']

        # Scaling factor T
        T = rl * (rg['count'] / gl['gene_length']).sum()

        # TPM
        if T == 0:
            tpm = pd.Series(0, index=rg.index)
        else:
            tpm = ((1e6 * rl) / float(T)) * (rg['count'] / gl['gene_length'])

        # Add sample column
        df[sample_name] = tpm

    # Output
    df.to_csv(sys.stdout, sep='\t')
    logging.info("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-n', '--sample_names', nargs='*',
                        help="Sample names, in the same order as coverage_files")
    parser.add_argument('-c', '--coverage_files', nargs='*',
                        help="Coverage files with 'sequence_id  count'")
    parser.add_argument('-i', '--sample_info',
                        help="TSV: sample_id   avg_read_length")
    parser.add_argument('-l', '--gene_lengths',
                        help="TSV: gene_id   gene_length")

    args = parser.parse_args()
    main(args)