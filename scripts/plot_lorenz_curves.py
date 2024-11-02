#!/usr/bin/env python3
'''
Plot Lorenz curves for a list of samples. Lorenz curves were first developed as a measure of income inequality. 
The dotted line represents perfectly uniform coverage. The x-intercept represents the fraction of the genome that 
lacks any coverage. Note coverage uniformity can only be compared fairly between samples with equal sequencing depth.

Input files were generated with Picard as follows, with a non-default value for COVERAGE_CAP of 500:

java -Xmx4g -jar picard.jar CollectWgsMetrics \
       -I <input.bam> \
       -O <output.txt> \
       -R <reference.fa> \
       -Q <min_base_quality> \
       -MQ <min_mapping_quality> \
       -CAP <coverage_cap>

Example usage:

python /projects/steiflab/research/hmacdonald/plots/plot_lorenz_curves.py \
       --sample_ids A138856_Iso5_Livecells SiHa_Livecells \
       --wgs_metrics_files Iso5_merged.wgs_metrics.txt SiHa_LiveCells.wgs_metrics.txt \
       --out_metrics SiHa_comparison_lorenz_out_metrics.csv \
       --out_hist SiHa_comparison_lorenz_out_histogram.csv \
       --out_plot SiHa_comparison_lorenz_out_plot.pdf

'''
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--sample_ids', nargs='+',
                        help='''List of sample identifiers corresponding to the WGS metrics files.''')

    parser.add_argument('--wgs_metrics_files', nargs='+',
                        help='''Path to one or more Picard WGS metrics files.''')

    parser.add_argument('--out_metrics',
                        help='''Path to output summary metrics file for all samples in .csv format.''')

    parser.add_argument('--out_hist',
                        help='''Path to output summary coverage histogram file for all samples in .csv format.''')

    parser.add_argument('--out_plot',
                        help='''Path to output Lorenz curve plot in .pdf format.''')

    args = parser.parse_args()

    return args




def read_wgs_metrics_file(sample_id, wgs_metrics_file):
    '''Parse the Picard WGS metrics file and format the genome metrics and coverage histogram.'''
    df = pd.read_table(os.path.join(wgs_metrics_file), sep='\t', comment='#')

    # create a separate data frame for the genome metrics and add the sample identifier
    df_metrics = df.iloc[[0]]
    df_metrics.insert(0, 'sample_id', sample_id)
    df_metrics.columns = [x.lower() for x in df_metrics.columns]

    # create a separate data frame for the coverage histogram and add the sample identifier
    # exclude the last row which represents loci with coverage greater than or equal to the coverage cap
    # Note: if the value in the last row is large, best to increase the coverage cap parameter 
    # when running Picard WGS metrics
    df_hist = df.iloc[2:,0:2]
    df_hist = df_hist.iloc[0:df_hist.shape[0]-1]
    df_hist.columns = ['coverage', 'count']
    df_hist = df_hist.set_index('coverage')
    df_hist.index.name = None
    df_hist = df_hist.T
    df_hist.index.name = 'sample_id'
    df_hist = df_hist.reset_index()
    df_hist['sample_id'][0] = sample_id

    return df_metrics, df_hist


def format_df_hist_row(df_hist_row):
    '''Reformat the histogram for a single sample.'''
    df_hist_row = pd.DataFrame(df_hist_row).T
    sample_id = df_hist_row.iloc[0]
    df_hist_row = df_hist_row.drop(['sample_id'], axis = 1)
    df_hist_row = df_hist_row.T
    df_hist_row = df_hist_row.reset_index()
    df_hist_row.columns = ['coverage', 'count']
    df_hist_row = df_hist_row.astype(int)

    return df_hist_row


def compute_cumulative_coverage(df_hist_row):
    '''Compute the cumulative coverage for a single sample for plotting.'''
    df_cumulative = format_df_hist_row(df_hist_row)
    print(df_hist_row)

    total_count = sum(df_cumulative['count'])
    df_cumulative['fraction_genome'] = df_cumulative['count'] / total_count
    df_cumulative['cumulative_fraction_genome'] = np.cumsum(df_cumulative['fraction_genome'])

    df_cumulative['total_coverage'] = df_cumulative['coverage'] * df_cumulative['count']
    total_coverage = sum(df_cumulative['total_coverage'])
    df_cumulative['fraction_total_coverage'] = df_cumulative['total_coverage'] / total_coverage
    df_cumulative['cumulative_fraction_total_coverage'] = np.cumsum(df_cumulative['fraction_total_coverage'])

    return df_cumulative



def plot_lorenz_curves(sample_ids, df_metrics, df_hist, out_plot, colour_list=None):
    '''Plot Lorenz curve for each sample, with coverage depth added to the label.
    sampleids = list of cells to plot '''
    if colour_list:
        cols = colour_list[0:len(sample_ids)]

    else:
        tableau_20 = ['#4e79a7', '#a0cbe8', '#f28e2b', '#ffbe7d', '#59a14f',
                      '#8cd17d', '#b6992d', '#f1ce63', '#499894', '#86bcb6',
                      '#e15759', '#ff9d9a', '#79706e', '#bab0ac', '#d37295',
                      '#fabfd2', '#b07aa1', '#d4a6c8', '#9d7660', '#d7b5a6']

        cols = tableau_20[0:len(sample_ids)]
    sns.set(context='talk',
            style='ticks',
            font='Helvetica',
            rc={'figure.titlesize': 15,
                'axes.titlesize': 15,
                'axes.labelsize': 15,
                'xtick.labelsize': 15,
                'ytick.labelsize': 15,
                'legend.fontsize': 15,
                'axes.linewidth': 1.3,
                'xtick.major.width': 1.3,
                'ytick.major.width': 1.3})

    fig = plt.figure(figsize=(6, 6))

    ax = fig.gca()

    for sample_id, col in zip(sample_ids, cols):
            # Check if sample_id is 'CellenONE' to change line style
        if sample_id == 'CellenONE':
            line_style = '--'  # Dotted line
        else:
            line_style = '-'   # Default solid line
        coverage_depth = df_metrics.loc[df_metrics['sample_id'] == sample_id, 'mean_coverage'].astype('float').item()

        sample_label = '{} ({:.2f}X)'.format(sample_id, coverage_depth)

        df_cumulative = compute_cumulative_coverage(df_hist[df_hist['sample_id'] == sample_id].transpose())

        ax.plot(df_cumulative['cumulative_fraction_genome'], 
                df_cumulative['cumulative_fraction_total_coverage'], 
                color=col, label=sample_label, linestyle=line_style)

    # add diagonal line representing perfectly uniform coverage
    plt.plot([0,1], [0,1], '--', color='0.1', linewidth=1)

    plt.xlim(0,1)
    plt.ylim(0,1)

    ax.set_xlabel('Cumulative fraction of genome')
    ax.set_ylabel('Cumulative fraction of total coverage')

    plt.legend(loc='upper left')

    sns.despine()
    plt.tight_layout()

    plt.savefig(out_plot)
    plt.close()



def main():
    args = parse_args()

    df_metrics = pd.DataFrame()
    df_hist = pd.DataFrame()

    for sample_id, wgs_metrics_file in zip(args.sample_ids, args.wgs_metrics_files):
        df_sample_metrics, df_sample_hist = read_wgs_metrics_file(sample_id, wgs_metrics_file)
        df_metrics = pd.concat([df_metrics, df_sample_metrics], ignore_index=True)
        df_hist = pd.concat([df_hist, df_sample_hist], ignore_index=True)

    df_metrics.to_csv(args.out_metrics, index=False)
    df_hist.to_csv(args.out_hist, index=False)

    cols = ['#a0cbe8', '#1f77b4', '#9467bd']
    plot_lorenz_curves(args.sample_ids, df_metrics, df_hist, args.out_plot, colour_list = cols)


if __name__ == '__main__':
    main()
