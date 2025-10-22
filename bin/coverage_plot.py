#!/usr/bin/env python3

"""
Generate coverage plots from depth files
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate coverage plots from depth files'
    )
    parser.add_argument(
        '-f', '--file',
        required=True,
        help='Input depth file (TSV format with columns: chrom, position, depth)'
    )
    parser.add_argument(
        '-t', '--title',
        default='Coverage Plot',
        help='Title for the plot'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output PDF file path'
    )
    parser.add_argument(
        '--xlabel',
        default='Position',
        help='X-axis label (default: Position)'
    )
    parser.add_argument(
        '--ylabel',
        default='Coverage',
        help='Y-axis label (default: Coverage)'
    )
    parser.add_argument(
        '--height',
        type=float,
        default=6.0,
        help='Figure height in inches (default: 6.0)'
    )
    parser.add_argument(
        '--aspect',
        type=float,
        default=2.0,
        help='Figure aspect ratio (default: 2.0)'
    )
    parser.add_argument(
        '--fill-color',
        default='skyblue',
        help='Fill color for area under curve (default: skyblue)'
    )
    parser.add_argument(
        '--line-color',
        default='navy',
        help='Line color (default: navy)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.5,
        help='Transparency of fill (0-1, default: 0.5)'
    )

    return parser.parse_args()


def read_depth_file(filepath):
    """Read depth file and return DataFrame"""
    try:
        df = pd.read_csv(filepath, sep='\t', header=None)
        df.columns = ["chrom", "position", "depth"]
        return df
    except Exception as e:
        sys.stderr.write(f"Error reading file {filepath}: {e}\n")
        sys.exit(1)


def create_coverage_plot(df, title, xlabel, ylabel, height, aspect,
                         fill_color, line_color, alpha, output):
    """Create and save coverage plot"""
    try:
        # Set seaborn theme
        sns.set_theme(style="darkgrid")

        # Create the plot
        fig = sns.relplot(
            x=df["position"],
            y=df["depth"],
            kind='line',
            height=height,
            aspect=aspect
        )

        # Add fill under the curve
        plt.fill_between(
            df["position"],
            df["depth"],
            0,
            facecolor=fill_color,
            color=line_color,
            alpha=alpha
        )

        # Set labels and title
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Save the figure
        plt.savefig(output, bbox_inches='tight')
        plt.close()

    except Exception as e:
        sys.stderr.write(f"Error creating plot: {e}\n")
        sys.exit(1)


def main():
    """Main function"""
    args = parse_args()

    # Read depth file
    df = read_depth_file(args.file)

    # Create and save plot
    create_coverage_plot(
        df=df,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        height=args.height,
        aspect=args.aspect,
        fill_color=args.fill_color,
        line_color=args.line_color,
        alpha=args.alpha,
        output=args.output
    )

    print(f"Coverage plot saved to: {args.output}")


if __name__ == "__main__":
    main()
