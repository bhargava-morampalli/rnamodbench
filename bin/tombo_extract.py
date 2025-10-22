#!/usr/bin/env python

"""
Extract raw statistics from Tombo stats files
Based on the Tombo Python API to extract regional statistics
"""

import argparse
import sys
import pandas as pd
from tombo import tombo_stats


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract statistics from Tombo stats files'
    )
    parser.add_argument(
        '--statfile',
        required=True,
        help='Input Tombo statistics file (.tombo.stats)'
    )
    parser.add_argument(
        '--chrom',
        required=True,
        help='Chromosome/contig name (e.g., 16s_88_rrsE or 23s_78_rrlB)'
    )
    parser.add_argument(
        '--strand',
        default='+',
        help='Strand (+ or -). Default: +'
    )
    parser.add_argument(
        '--start',
        type=int,
        required=True,
        help='Start position (1-based)'
    )
    parser.add_argument(
        '--end',
        type=int,
        required=True,
        help='End position (1-based, inclusive)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output CSV file path'
    )

    return parser.parse_args()


def extract_tombo_stats(statfile, chrom, strand, start, end, output):
    """
    Extract regional statistics from Tombo stats file
    
    Args:
        statfile: Path to Tombo statistics file
        chrom: Chromosome/contig name
        strand: Strand ('+' or '-')
        start: Start position (1-based)
        end: End position (1-based, inclusive)
        output: Output CSV file path
    """
    try:
        # Load Tombo statistics
        sample_level_stats = tombo_stats.LevelStats(statfile)
        
        # Get regional statistics
        reg_level_stats = sample_level_stats.get_reg_stats(chrom, strand, start, end)
        
        # Convert to DataFrame and save as CSV
        df = pd.DataFrame(reg_level_stats)
        df.to_csv(output, index=False)
        
        print(f"Successfully extracted {len(df)} positions to {output}")
        
    except Exception as e:
        sys.stderr.write(f"Error extracting Tombo statistics: {e}\n")
        sys.exit(1)


def main():
    """Main function"""
    args = parse_args()
    
    extract_tombo_stats(
        statfile=args.statfile,
        chrom=args.chrom,
        strand=args.strand,
        start=args.start,
        end=args.end,
        output=args.output
    )


if __name__ == "__main__":
    main()
