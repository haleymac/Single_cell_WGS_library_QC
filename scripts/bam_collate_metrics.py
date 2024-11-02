'''Collate and summarize metrics across cells.'''
import numpy as np
import pandas as pd
import os
from dictionary_functions import *


def _make_fastq_dictionaries(input_dir):
    '''make dictionaries (for paired end reads) with "cellids: path to fastq_screen output file
    input_dir: directory with fastq_screen output'''
    fastq_screen_dict_1 = {}
    fastq_screen_dict_2 = {}
    for file in os.listdir(input_dir):
        if file.endswith('screen.txt'):
            fq = file.split('_trimmed')[0]
            read_number = fq.split('_')[-1]
            print(read_number)
            cell_id = fq.split('_')[0]
            if read_number == '1':
                fastq_screen_dict_1[cell_id] = os.path.join(input_dir, file)
            elif read_number == '2':
                fastq_screen_dict_2[cell_id] = os.path.join(input_dir, file)
    return fastq_screen_dict_1, fastq_screen_dict_2



def parse_fastq_screen_file(cell_id, fastq_screen_file, read=None):
    '''Parse a fastq screen txt file for a single FASTQ file.'''
    df = pd.read_csv(fastq_screen_file, sep='\t', comment='#', header=None, skiprows=2)

    columns = ['genome', 'reads_processed', 'unmapped', 'percent_unmapped',
            'one_hit_one_genome', 'percent_one_hit_one_genome',
            'multiple_hits_one_genome', 'percent_multiple_hits_one_genome',
            'one_hit_multiple_genomes', 'percent_one_hit_multiple_genomes', 
            'multiple_hits_multiple_genomes', 'percent_multiple_hits_multiple_genomes']

    df.columns = columns

    hit_no_genomes_bool = df['genome'].str.contains('%Hit_no_genomes:')

    percent_hit_no_genomes = float(df[hit_no_genomes_bool]['genome'].iloc[0].split(': ')[1])

    df_hit_no_genomes = pd.DataFrame({'genome': 'No hit',
                                    'variable': 'percent_hit_no_genomes',
                                    'value': percent_hit_no_genomes}, index=[0])

    df = df.drop(df[hit_no_genomes_bool].index, axis=0)

    df_melt = df.melt(id_vars='genome')

    df_melt = pd.concat([df_melt, df_hit_no_genomes], axis=0)

    if read is not None:
        df_melt.insert(0, 'read', read)

    df_melt.insert(0, 'cell_id', cell_id)

    return df_melt



def collate_fastq_screen_single_end(fastq_screen_dict, out_file):
    '''Collate the fastq screen output from multiple single-end FASTQ files into a single melted data frame.'''
    df_collate = pd.DataFrame()

    for cell_id in fastq_screen_dict.keys():
        df_melt = parse_fastq_screen_file(cell_id, fastq_screen_dict[cell_id], read='R1')
        df_collate = pd.concat([df_collate, df_melt], axis=0)

    df_collate.to_csv(out_file, index=False)


def collate_fastq_screen_paired_end(input_dir, out_file):
    '''Collate the fastq screen output from multiple paired-end FASTQ files into a single melted data frame.'''
    
    fastq_screen_dict_1, fastq_screen_dict_2 = _make_fastq_dictionaries(input_dir) 
    print(fastq_screen_dict_1)

    df_collate = pd.DataFrame()

    for cell_id in fastq_screen_dict_1.keys():
        if cell_id in fastq_screen_dict_2.keys():
            try:
                df_melt_1 = parse_fastq_screen_file(cell_id, fastq_screen_dict_1[cell_id], read='R1')
            except pandas.errors.EmptyDataError:
                print(cell_id + " is empty")
            try:
                df_melt_2 = parse_fastq_screen_file(cell_id, fastq_screen_dict_2[cell_id], read='R2')
            except pandas.errors.EmptyDataError:
                print(cell_id + " is empty")
            df_collate = pd.concat([df_collate, df_melt_1, df_melt_2], axis=0)
    df_collate.to_csv(out_file, index=False)


indir = '/projects/steiflab/research/hmacdonald/isolatrix/isolatrix_scratch/data_scratch/IX12003_scratch/fastq_metrics_output/fastq_screen_output/individual_output'
out = '/projects/steiflab/research/hmacdonald/isolatrix/isolatrix_scratch/data_scratch/IX12003_scratch/fastq_metrics_output/fastq_screen_output/test.csv'
#collate_fastq_screen_paired_end(indir, out)



def parse_duplicate_metrics(cell_id, duplicate_metrics_file):
    '''Parse the Picard duplicate metrics file.'''
    df = pd.read_csv(duplicate_metrics_file, sep='\t', comment='#', nrows=1)

    df.columns = [x.lower() for x in df.columns.values]

    df = df.rename(columns={'percent_duplication': 'duplicate_fraction'})

    df.insert(0, 'cell_id', cell_id)

    return df


def parse_mt_coverage_bed(cell_id, mt_coverage_file):
    df = pd.read_csv(mt_coverage_file, sep = '\t', names = ['chr', 'start', 'end', 'mt_coverage_depth', 'n_bases_covered', 'len_coverage', 'fraction_covered'])
    df.insert(0, 'cell_id', cell_id)
    df = df[['cell_id', 'chr', 'start', 'end', 'mt_coverage_depth']]
    return df


def collate_mt_coverage(mt_coverage_dict, outfile):
    df_collate = pd.DataFrame()
    for cell_id in mt_coverage_dict.keys():
        try:
            df = parse_mt_coverage_bed(cell_id, mt_coverage_dict[cell_id])
            df_collate = pd.concat([df_collate, df], axis=0)
        except pd.errors.EmptyDataError:
            print(f"{cell_id} has an empty mt coverage file")
    df_collate.to_csv(outfile, index = False)
    
    
def collate_duplicate_metrics(duplicate_metrics_dict, out_file):
    '''Collate the Picard duplicate metrics output from multiple cells.'''
    df_collate = pd.DataFrame()

    for cell_id in duplicate_metrics_dict.keys():
        df = parse_duplicate_metrics(cell_id, duplicate_metrics_dict[cell_id])
        df_collate = pd.concat([df_collate, df], axis=0)
            
    df_collate.to_csv(out_file, index=False)



def compute_coverage_breadth(df_hist_col):
    '''Compute the coverage breadth from the WGS metrics historgram.'''
    genome_territory = sum(df_hist_col['count'])
    empty_territory = df_hist_col['count'][0]

    covered_count = genome_territory - empty_territory
    coverage_breadth = covered_count / genome_territory

    return coverage_breadth


def parse_wgs_metrics(cell_id, wgs_metrics_file):
    '''Parse the Picard WGS metrics file.'''
    df = pd.read_csv(wgs_metrics_file, sep='\t', comment='#')

    df_metrics = df.iloc[[0]]
    df_metrics.insert(0, 'cell_id', cell_id)
    df_metrics.columns = [x.lower() for x in df_metrics.columns]

    df_hist_col = df.iloc[2:, 0:2].astype(int).reset_index(drop=True)
    df_hist_col.columns = ['coverage', 'count']

    coverage_breadth = compute_coverage_breadth(df_hist_col)
    df_metrics.insert(2, 'coverage_breadth', coverage_breadth)

    df_hist = df_hist_col.set_index('coverage')
    df_hist.index.name = None
    df_hist = df_hist.T
    df_hist = df_hist.reset_index(drop=True)
    df_hist.insert(0, 'cell_id', cell_id)

    return df_metrics, df_hist



def parse_rna_metrics(cell_id, rna_metrics_file):
    df = pd.read_csv(rna_metrics_file, sep='\t', comment='#')
    df_metrics = df.iloc[[0]]
    df_metrics.insert(0, 'cell_id', cell_id)
    df_metrics.columns = [x.lower() for x in df_metrics.columns]
    
    
    df_hist_col = df.iloc[2:, 0:2].astype(float).reset_index(drop=True)
    df_hist_col.columns = ['normalized_position', 'normalized_coverage']
    df_hist = df_hist_col.set_index('normalized_position')
    df_hist.index.name = None
    df_hist = df_hist.T
    df_hist = df_hist.reset_index(drop=True)
    df_hist.insert(0, 'cell_id', cell_id)
    
    return df_metrics, df_hist




def collate_rna_metrics(rna_metrics_dict, out_file, out_hist):
    '''Collate the Picard RNA metrics output from multiple cells.'''
    df_metrics_collate = pd.DataFrame()
    df_hist_collate = pd.DataFrame()
    
    for cell_id in rna_metrics_dict.keys():
        df_metrics, df_hist = parse_rna_metrics(cell_id, rna_metrics_dict[cell_id])
        df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
        df_hist_collate = pd.concat([df_hist_collate, df_hist], axis=0)
    
    df_metrics_collate.to_csv(out_file, index=False)
    df_hist_collate.to_csv(out_hist, index=False)


def collate_wgs_metrics(wgs_metrics_dict, out_file, out_hist):
    '''Collate the Picard WGS metrics output from multiple cells.'''
    df_metrics_collate = pd.DataFrame()
    df_hist_collate = pd.DataFrame()

    for cell_id in wgs_metrics_dict.keys():
        df_metrics, df_hist = parse_wgs_metrics(cell_id, wgs_metrics_dict[cell_id])
        df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
        df_hist_collate = pd.concat([df_hist_collate, df_hist], axis=0)

    df_metrics_collate.to_csv(out_file, index=False)
    df_hist_collate.to_csv(out_hist, index=False)


def parse_gc_bias_metrics(cell_id, gc_bias_metrics_file):
    '''Parse the Picard GC bias metrics file.'''
    df = pd.read_csv(gc_bias_metrics_file, sep='\t', comment='#')

    df.columns = [x.lower() for x in df.columns]

    df.insert(0, 'cell_id', cell_id)

    return df


def collate_gc_bias_metrics(gc_bias_metrics_dict, gc_bias_summary_dict, out_file, out_summary):
    '''Collate the Picard GC bias metrics output from multiple cells.'''
    df_metrics_collate = pd.DataFrame()
    df_summary_collate = pd.DataFrame()
    emptycellsFile = open('cells_with_empty_gc_metrics_files.txt', 'w')
    print('The follow cells have not been included in the collated gc metrics csv, because their gc metrics txt files were empty', file = emptycellsFile)

    for cell_id in gc_bias_metrics_dict.keys():
        try:
            df = pd.read_csv(gc_bias_metrics_dict[cell_id], sep='\t', comment='#')
            df_metrics = parse_gc_bias_metrics(cell_id, gc_bias_metrics_dict[cell_id])
            df_summary = parse_gc_bias_metrics(cell_id, gc_bias_summary_dict[cell_id])
            df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
            df_summary_collate = pd.concat([df_summary_collate, df_summary], axis=0)

        except pd.errors.EmptyDataError:
            print(cell_id, file = emptycellsFile)
            
    df_metrics_collate.to_csv(out_file, index=False)
    df_summary_collate.to_csv(out_summary, index=False)

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def parse_insert_metrics_old(cell_id, insert_metrics_file):
    '''Parse the Picard duplicate metrics file.'''
    df = pd.read_csv(insert_metrics_file, sep='\t', comment='#')
    print('here')
    df_metrics = df.iloc[[0]]
    df_metrics.insert(0, 'cell_id', cell_id)
    df_metrics.columns = [x.lower() for x in df_metrics.columns]
    df_tandem = df.iloc[[0, 1, 2]]
    #this checks if you have tandem and RF reads in sample (if row 4 of df had tandem or rf reads)
    if df_tandem.iloc[2, 8] == 'TANDEM' or df_tandem.iloc[2, 8] == 'RF':
        row1_fraction = [df_tandem.iloc[1,7]/(df_tandem.iloc[0,7] + df_tandem.iloc[2,7])]
        row2_fraction = [df_tandem.iloc[2,7]/(df_tandem.iloc[0,7] + df_tandem.iloc[1,7])]
        if df_tandem.iloc[1,8] == 'TANDEM':
            df_metrics['tandem_fraction'] = row1_fraction
        elif df_tandem.iloc[1,8] == 'RF':
            df_metrics['rf_fraction'] = row1_fraction
            
        if df_tandem.iloc[2,8] == 'TANDEM':
            df_metrics['tandem_fraction'] = row2_fraction
        elif df_tandem.iloc[2,8] == 'RF':
            df_metrics['rf_fraction'] = row2_fraction
    #else if you only have tandem or rf reads and not both
    else: 
        tandem_rf_fraction = [df_tandem.iloc[1,7]/df_tandem.iloc[0,7]]
        if df_tandem.iloc[1,8] == 'TANDEM':
            df_metrics['tandem_fraction'] = tandem_rf_fraction
            df_metrics['rf_fraction'] = [0]
        elif df_tandem.iloc[1,8] == 'RF':
            df_metrics['rf_fraction'] = tandem_rf_fraction
            df_metrics['tandem_fraction'] = [0]
        else:
            df_metrics['tandem_fraction'] = [0]
            df_metrics['rf_fraction'] = [0]

    if df_metrics.tandem_fraction[0] == 0 and df_metrics.rf_fraction[0] == 0:
        hist_start = 1
    elif (df_metrics.tandem_fraction[0] == 0 and df_metrics.rf_fraction[0] != 0) or (df_metrics.tandem_fraction[0] != 0 and df_metrics.rf_fraction[0] == 0):
        hist_start = 2
    else:
        hist_start = 3
    
    df_hist_col = df.iloc[hist_start+1:, 0:2].astype(int).reset_index(drop=True)
    df_hist_col.columns = ['insert_size', 'read_count']
    df_hist = df_hist_col.set_index('insert_size')
    df_hist.index.name = None
    df_hist = df_hist.T
    df_hist = df_hist.reset_index(drop=True)
    df_hist.insert(0, 'cell_id', cell_id)
    return df_metrics, df_hist


def _find_section_start_end_lines(file_path, section_header):
    lines = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith(section_header):
                lines.append(i)
            elif (len(lines) == 1) and (line.strip() == ''):
                lines.append(i) 
    return lines


def parse_insert_metrics(cell_id, insert_metrics_file):
    '''Parse the Picard duplicate metrics file.'''
    metrics_positions = _find_section_start_end_lines(insert_metrics_file, 'MEDIAN_INSERT_SIZE')
    metric_nrows = (metrics_positions[1]) - (metrics_positions[0])
    metrics_df = pd.read_csv(insert_metrics_file, sep='\t', skiprows=metrics_positions[0], nrows=metric_nrows-1)
    metrics_df.insert(0, 'cell_id', cell_id)
    metrics_df.columns = [x.lower() for x in metrics_df.columns]
    
    histogram_positions = _find_section_start_end_lines(insert_metrics_file, 'insert_size')
    hist_nrows = histogram_positions[1] - histogram_positions[0]
    hist_df = pd.read_csv(insert_metrics_file, sep='\t', skiprows=histogram_positions[0], nrows=hist_nrows)
    hist_df = hist_df.astype(int).reset_index(drop=True)
    hist_df.columns = ['insert_size', 'read_count']
    hist_df = hist_df.set_index('insert_size')
    hist_df.index.name = None
    hist_df = hist_df.T
    hist_df = hist_df.reset_index(drop=True)
    hist_df.insert(0, 'cell_id', cell_id)
    
    return metrics_df, hist_df


def collate_insert_metrics(insert_metrics_dict, out_file, out_hist):
    df_metrics_collate = pd.DataFrame()
    df_hist_collate = pd.DataFrame()
    for cell_id in insert_metrics_dict.keys():
        if is_non_zero_file(insert_metrics_dict[cell_id]) == True:
            #insert_df = pd.read_csv(insert_metrics_dict[cell_id], sep='\t', comment='#')
            df_metrics, df_hist = parse_insert_metrics(cell_id, insert_metrics_dict[cell_id])
            df_metrics_collate = pd.concat([df_metrics_collate, df_metrics], axis=0)
            df_hist_collate = pd.concat([df_hist_collate, df_hist], axis=0)
        df_metrics_collate.to_csv(out_file, index=False)
        df_hist_collate.to_csv(out_hist, index=False)



"""
cells = []
for file in os.listdir('/projects/steiflab/research/hmacdonald/total_RNA/data/buettner/buettner_scratch/bams/star_basic_alignment'):
    if file.endswith('.bam'):
        cell = file.split('.')[0]
        if cell not in cells:
            cells.append(cell)
    
output_dict = create_output_file_dict(cells, '/projects/steiflab/research/hmacdonald/total_RNA/data/buettner/buettner_scratch/bams/bam_metrics_output_no_dups')

    
metrics_df, hist_df = parse_insert_metrics('ERR580236', output_dict['output']['insert']['ERR580236'])
print(hist_df)
"""