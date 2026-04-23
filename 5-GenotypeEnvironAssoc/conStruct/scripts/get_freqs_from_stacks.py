#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract allele frequencies from Stacks sumstats")
    parser.add_argument('--sumstats', required=True, help='Path to Stacks .p.sumstats.tsv file')
    parser.add_argument('--out', required=True, help='Path to output frequencies CSV file')
    args = parser.parse_args()

    # Read the sumstats file, skipping comment lines starting with '#'
    stacks_res = pd.read_csv(args.sumstats, sep='\t', comment='#')
    
    # In older stacks versions or depending on header structure, Locus ID column might be 'Locus ID' or '# Locus ID'
    locus_col = '# Locus ID' if '# Locus ID' in stacks_res.columns else 'Locus ID'
    if locus_col not in stacks_res.columns and 'Locus ID' in stacks_res.columns:
        locus_col = 'Locus ID'

    # Pivot to get frequencies
    freqs = stacks_res[[locus_col, 'Pop ID', 'P']].pivot(index='Pop ID', columns=locus_col, values='P')
    
    # Save to csv. index_label=False prevents adding a header for the index column, 
    # making it compatible with R's read.csv(..., row.names=1)
    freqs.to_csv(args.out, index_label=False)

if __name__ == '__main__':
    main()
