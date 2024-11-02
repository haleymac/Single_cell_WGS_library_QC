import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.lines as mlines
import argparse
import seaborn as sns
import sys


'''
This script contains functions to make plots representing different QC metrics. It has individual functions that can be called to output one .png plot.
It has other functions that produce particular helpful plots for particular collated bam_metrics output files (from bam_metrics pipeline). 
It had one large function to make many detailed library plots for all collated bam_metrics output files as seperate pdf files.
It has another function that outputs basic standard plots to get an overview of library quality in one pdf
'''

def parse_args():
    '''Parse command line arguments '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir',
                        help = '''Directory to save pdf files containing plots in.''')
    parser.add_argument('--collated_input_dir',
                        help = '''Directory containg collated gatk metric files output by the bam_metrics pipeline''')
    parser.add_argument('--iso_csv',
                        help = '''Path to the isolatrix/sequencing csv provided by Eric''')
    parser.add_argument('--collated_fastq_screen_csv',
                        help = '''Path to the collated fastq screen metrics csv output by the fastq_metrics pipeline''')

    args = parser.parse_args()
    return args


#set colour palette and theme for plots
pal = sns.color_palette()
sns.set_palette(pal)
sns.set_context("talk")

#set standard colours for particular cell conditions - note OddCell gets collapsed into livecell condition currently 
colour_dict = {'LiveCell': '#1f77b4',
                'NCC': '#d62728',
                'NTC': '#ff7f0e',
                'gDNA': '#2ca02c', 
                'OddCell': '#9467bd'}

seaborn_colours = {'blue': '#1f77b4', 
                   'orange': '#ff7f0e', 
                   'green': '#2ca02c', 
                   'red' : '#d62728', 
                   'purple': '#9467bd', 
                   'brown': '#8c564b', 
                   'pink': '#e377c2',
                   'gray': '#7f7f7f', 
                   'yellow': '#bcbd22', 
                   'light_blue': '#17becf'}

def make_colour_dict(value_lst): 
    '''make a dictionary with values assigned to the default seaborn palette - useful for functions when you want to change what a plot is colouring by 
    value_lst: should be the uniue values of a df column that you want to colour a plot by'''
    pal = sns.color_palette().as_hex()
    i = 0
    colour_dict = {}
    for key in value_lst:
        colour_dict[key] = pal[i]
        i+=1
    return colour_dict
    
def make_metric_name_look_nice(metric):
    '''changes the name of the metric you want to plot - removes underscores and capitolizes the first letter - use this when labelling axis'''
    dont_change_list = ['gDNA', 'NCC', 'NTC']
    split_metric = metric.split('_')
    nicer_metric_name = ''
    for i, word in enumerate(split_metric):
        if i == 0 and word not in dont_change_list:
            metric_word = word.capitalize()
        else:
            metric_word = word
        nicer_metric_name = str(nicer_metric_name + ' ' + metric_word)
    return nicer_metric_name    
    
def make_palette(df, condition):
    '''make a new palette using seaborn default colours
    df: df to make the palette for
    condition: df column with conditions to make palette for (ex. cell_condition with 'LiveCell' 'gDNA' etc)'''
    new_pal = {}
    for cond, colour in zip(df[condition].unique(), seaborn_colours.keys()):
        new_pal[cond] = seaborn_colours[colour]
    return new_pal


def merge_sequencing_csvs_with_metrics_csvs(collated_input_dir, sequencing_csv, chip_id, iso_library = True, paired_end = True):
    '''reads in the csv file that comes with the sequencing run and combines row, col, library and pool columns to create a cellid that matches the output from the bam_metrics pipeline.
    Subsequently merges the new sequencing df with the collated metrics csvs produced by the bam metrics pipeline
    Note: doesn't merge the histogram csvs, only the summary metrics csvs
    Collated_input_dir: is the path to the directory with the collated bam metrics output
    sequencing_df: the dataframe produced by running the sequencing 
    Returns a dictionary where each key is a particular print chip, with the values being the df produced by combining the collated metrics from the bam_metrics pipeline.
    '''
    #metrics from collated bam metrics files to include in the new df
    desired_metrics = {'duplicate_metrics':['cell_id','unpaired_reads_examined', 'read_pairs_examined','unmapped_reads','duplicate_fraction','estimated_library_size', 'unpaired_read_duplicates', 'read_pair_duplicates'], 
                    'gc_bias_summary':['cell_id','total_clusters', 'aligned_reads', 'at_dropout', 'gc_dropout', 'gc_nc_0_19', 'gc_nc_20_39', 'gc_nc_40_59',
       'gc_nc_60_79', 'gc_nc_80_100'], 
                    'insert_metrics':['cell_id','median_insert_size','min_insert_size','max_insert_size','mean_insert_size','standard_deviation', 'tandem_fraction', 'rf_fraction'],
                    'wgs_metrics':['cell_id','coverage_breadth','mean_coverage','sd_coverage','median_coverage', 'pct_1x', 'pct_5x', 'pct_10x', 'pct_15x', 'pct_20x', 'pct_25x', 'pct_30x', 'pct_40x', 'pct_50x', 'pct_60x', 'pct_70x', 'pct_80x', 'pct_90x', 'pct_100x'],
                    'mt_coverage': ['cell_id', 'mt_coverage_depth'],
                    'RNA_Metrics': ['cell_id' ,'pf_bases','pf_aligned_bases', 'ignored_reads','correct_strand_reads','incorrect_strand_reads','num_r1_transcript_strand_reads','num_r2_transcript_strand_reads','num_unexplained_reads','pct_r1_transcript_strand_reads','pct_r2_transcript_strand_reads','pct_ribosomal_bases','pct_coding_bases','pct_utr_bases','pct_intronic_bases','pct_intergenic_bases','pct_mrna_bases','pct_usable_bases','pct_correct_strand_reads','median_cv_coverage','median_5prime_bias','median_3prime_bias','median_5prime_to_3prime_bias','sample','library','read_group', 'coding_bases', 'utr_bases', 'intronic_bases']}
    iso_df = pd.read_csv(sequencing_csv) 
    
    if iso_library == True:
        #create cell_id column in isolatrix/sequencing df - can't do this for other libraries 
        iso_df['cell_id'] = (iso_df['chip_id'] + '-' + iso_df['pool_id'] + '-R' + iso_df['row'].apply(str) + '-C' + iso_df['column'].apply(str))

    iso_df_merged = iso_df
    present_metrics_files = []
    i = 0
    for file in os.listdir(collated_input_dir):
        chip = file.split('.')[0]
        metric_class = file.split('.')[1]
        if chip == chip_id:
            if metric_class in desired_metrics.keys():
                present_metrics_files.append(metric_class)
                df = pd.read_csv(os.path.join(collated_input_dir, file))
                df = df.loc[:,desired_metrics[metric_class]]
                iso_df_merged = iso_df_merged.merge(df, on='cell_id', how='left')
                i += 1
    if i == 0:
        print('No metrics files found - check chip id given to function (should be the chip id from the metrics files (i.e. chip_id.metricfile.csv))')
    else:
        for expected_metric_file in desired_metrics.keys():
            if expected_metric_file not in present_metrics_files:
                print(f'Missing collated metrics file for {expected_metric_file}')
        #calculate total raw reads and mapping percent from gatk insert metrics file
    iso_df_merged['n_unique_molecules'] = ((iso_df_merged['unpaired_reads_examined'].astype('float') - iso_df_merged['unpaired_read_duplicates'].astype('float')) + (iso_df_merged['read_pairs_examined'].astype('float') - iso_df_merged['read_pair_duplicates'].astype('float')))
    iso_df_merged['total_raw_reads'] = (iso_df_merged['unpaired_reads_examined'].astype('float') + iso_df_merged['read_pairs_examined'].astype('float')*2 + iso_df_merged['unmapped_reads'].astype('float'))
    iso_df_merged['mapping_percent'] = ((iso_df_merged['unpaired_reads_examined'].astype('float') + iso_df_merged['read_pairs_examined'].astype('float')*2)/iso_df_merged['total_raw_reads'])*100
    iso_df_merged['mapping_percent'] = (iso_df_merged['mapping_percent'].astype('float'))
    iso_df_merged = iso_df_merged.rename({'mean_coverage': 'mean_coverage_depth'}, axis=1)
    #older libraries don't have mt_coverage files so we use try/except here 
    try:
        iso_df_merged['mt_coverage_depth'] = (iso_df_merged['mt_coverage_depth'].astype('float'))
        iso_df_merged['mt_read_percent'] = (iso_df_merged['mt_coverage_depth']/(iso_df_merged['read_pairs_examined']*2))*100
    except KeyError:
        print('No file for mt coverage')
    #create new column to indicate whether a LiveCell initially was labelled as an OddCell 
    #iso_df_merged['odd_cell'] = np.where(iso_df_merged['cell_condition'] == 'OddCell', True, False)
    #collapse OddCell condition into 'LiveCell condition
    #iso_df_merged.replace('OddCell', 'LiveCell', inplace = True)
    #add print chip to the chip dict library for this sequencing run 
    return iso_df_merged


def make_run_metrics_table(collated_metrics_df, outdir, save_pdf = 'No'):
    '''Combines collated run metrics for each cell into one dataframe. Performs revelant operations (mean, median, max, sum) to summarize overall run metrics and outputs a nice table.
    collated_bam_metrics_output_dir: the directory with collated output from bam_metrics pipeline
    collated_metrics_df: df for a particular chip produced by merge_sequencing_csvs_with_metrics_csvs
    outdir: place to save resulting table
    save_pdf: if 'No' will save table as a png in outdir, if anything else, will save in a pdf file with the name of the argument'''
    
    collated_df = collated_metrics_df[['cell_condition', 'chip_id', 'total_raw_reads', 'read_pairs_examined',
    'unmapped_reads', 'coverage_breadth', 'mean_coverage',
        'mapping_percent', 'duplicate_fraction',
    'median_insert_size', 'min_insert_size', 'max_insert_size',
    'tandem_fraction', 'rf_fraction']]
    chip_id = collated_df.chip_id[0]
    df_summary = collated_df.groupby('cell_condition').agg({'total_raw_reads':['sum'], 
                                                            'read_pairs_examined' : ['sum'], 
                                                            'unmapped_reads' : ['sum'], 
                                                            'coverage_breadth' : ['mean'], 
                                                            'mean_coverage': ['mean'],
                                                            'mapping_percent': ['mean'],
                                                            'duplicate_fraction' : ['mean'],
                                                            'median_insert_size' : ['median'],
                                                            'min_insert_size': ['min'], 
                                                            'max_insert_size': ['max'],
                                                            'tandem_fraction': ['mean'], 
                                                            'rf_fraction': ['mean']})
    df_summary = df_summary.round(3)
    df_summary = df_summary.transpose()
    plt.rcParams['savefig.dpi'] = 300
    row_labels = ['Sum total raw reads', 'Sum read pairs',
        'Sum unmapped_reads', 'Mean coverage breadth', 'Mean coverage depth',
            'Mean mapping percent', 'Mean duplicate fraction',
        'Median insert size', 'Minimum insert size', 'Maximum insert size',
        'Mean tandem reads fraction', 'Mean RF reads fraction']
    rcolors = plt.cm.BuPu(np.full(len(row_labels), 0.1))
    ccolors = plt.cm.BuPu(np.full(len(df_summary.columns), 0.1))
    fig, ax = plt.subplots()
    # hide axes
    fig.patch.set_visible(False)
    ax.axis('off')

    table = ax.table(cellText=df_summary.values, 
                    colLabels=df_summary.columns, 
                    rowLabels = row_labels, 
                    rowColours=rcolors,
                    rowLoc='right',
                    colColours=ccolors,
                    loc='center')
    table.scale(0.5, 0.5)
    fig.tight_layout()
    plt.title(f'{chip_id} metrics summary')
            
    if save_pdf == 'No':
        plt.savefig(os.path.join(outdir, f'{chip_id}_metrics_table.png'), bbox_inches = 'tight')
        plt.close('all')
        
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')        
            


def make_violin(df, split_metric, plot_metric, title = 'default', legend = True, handles = 'default', bbox = (1,1)):
    '''make a violin plot to compare distribution of a metric between different groups. 
    df = dataframe to plot
    split_metric: df column of groups you want to compare - is x axis variable
    plot_metric: df column of the metric you want to compare - is y axis variable
    legend: set True to create a legend for different groups
    handles: pass custom handles for the legend
    bbox: coordinates to place legend at'''
    cut_at = df[plot_metric].min()
    categories = df[split_metric].unique()
    for category in categories:
        len_df_category = len(df[df[split_metric] == category])
        df.replace(category, f'{category} (n = {len_df_category})', inplace = True)
    plt.figure()
    vi = sns.violinplot(data=df, y = plot_metric, x = split_metric,  color = 'lightgrey', cut = cut_at)
    if legend == True:
        sw = sns.stripplot(data=df, y = plot_metric, x = split_metric, hue = split_metric, size = 3, edgecolor = 'black')
        if handles == 'default':
            sw.legend(bbox_to_anchor=bbox)
        else:
            sw.legend(handles = handles, bbox_to_anchor=bbox)
    else:
        sw = sns.stripplot(data=df, y=plot_metric, x = split_metric, hue = split_metric, size = 3, edgecolor = 'black', legend = False)
    plt.xlabel(' ')
    plt.ylabel(make_metric_name_look_nice(plot_metric))
    if title == 'default':
        plt.title(f'{make_metric_name_look_nice(plot_metric)}: {make_metric_name_look_nice(split_metric)}')
    else:
        plt.title(title)
    plt.savefig(f'{split_metric}_{plot_metric}_violin.png',  bbox_inches = 'tight', dpi = 600)


def make_swarm_boxplot(df, xmetric, ymetric, title = 'default', log = False, legend = True, include_n_in_legend = True, handles = 'default', bbox = (1,1), outplot = 'swarm_boxplot.png', palette = colour_dict, xlabel = ' ', ylabel = 'default', x_rotate = False, include_n = False, extra_code=None):
    #df[xmetric] = df[ymetric].
    categories = df[xmetric].unique()
    if include_n == True:
        for category in categories:
            len_df_category = len(df[df[xmetric] == category])
            df.replace(category, f'{category} (n = {len_df_category})', inplace = True)
    if include_n == True and include_n_in_legend == True:
        palette = update_pal_with_n(df, condition = xmetric, palette = palette)

    if include_n == False and include_n_in_legend == True:
        handles = make_handles(df, xmetric, palette = palette, include_n = True, include_median = False)

    fig, ax = plt.subplots()
    if log == True:
        df[ymetric] = df[ymetric] + 0.01
        ax.set_yscale('log')
        #ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x-1)))
        #ax.yaxis.set_major_locator(ticker.FixedLocator(ymetric))
    vi = sns.boxplot(data=df, y=ymetric, x = xmetric,  color = 'lightgrey')
    if legend == True:
        if palette == False:
            sw = sns.stripplot(data=df, x = xmetric, y = ymetric, hue = xmetric, size = 3, edgecolor = 'black')
        else:
            sw = sns.stripplot(data=df, y=ymetric, x = xmetric, hue = xmetric, palette = palette, size = 3, edgecolor = 'black')
        if handles == 'default':
            sw.legend(bbox_to_anchor=bbox)
        else:
            sw.legend(handles = handles, bbox_to_anchor=bbox)
    else:
        if palette == False:
            sw = sns.stripplot(data=df, y=ymetric, x = xmetric, hue = xmetric, size = 3, edgecolor = 'black', legend = False)
        else:
            sw = sns.stripplot(data=df, y=ymetric, x = xmetric, hue = xmetric, palette = palette, size = 3, edgecolor = 'black', legend = False)
    
    plt.xlabel(xlabel)
    if x_rotate == True:
        plt.xticks(rotation = 45)
    if ylabel == 'default':
        plt.ylabel(make_metric_name_look_nice(ymetric))
    else:
        plt.ylabel(ylabel)   
    if title == 'default':
        plt.title(f'{make_metric_name_look_nice(ymetric)}: {make_metric_name_look_nice(xmetric)}')
    else:
        plt.title(title)
    if extra_code:
        exec(extra_code, globals(), locals())
    plt.savefig(outplot,  bbox_inches = 'tight', dpi = 600)
    plt.close('all')
    
    
def _make_faceted_combined_boxplots(dataframe, metric, outdir, facet_by = 'pool_id', coloured_by = 'cell_condition', log = 'yes', save_pdf='No', bbox = (1,1)):
    '''makes one combined, faceted boxplot by a given category (ex. pool_id or chip_id). Boxplot is made based on given metric.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: QC metric you want boxplots for
    outdir: where you want the plot to be generated
    facet_by = category to facet by (ex. pool_id or chip_id)
    coloured_by: boxplots are automatically coloured by cell_condition, but you can specify something else if you want'''
    chipid_list = list(dataframe.chip_id.unique())
    for chip in chipid_list:
        plt.figure()
        ax = sns.boxplot(x=facet_by, y=metric, hue = 'cell_condition', palette = colour_dict, data=dataframe)
        plt.ylabel(make_metric_name_look_nice(metric))
        plt.xlabel(make_metric_name_look_nice(facet_by))
        plt.title(f"{chip} {make_metric_name_look_nice(metric)} for each {(facet_by).replace('_', ' ')}")
        plt.ticklabel_format(style='plain', axis='y')
        if log == 'yes':
            plt.yscale('log')
        df_num_cells = dataframe.drop_duplicates().groupby('cell_condition').size().reset_index() 
        handles = []
        for index, row in df_num_cells.iterrows():
            handle = mlines.Line2D([], [], color=colour_dict[row['cell_condition']],
                        markersize=25, label=str(f"{row['cell_condition']} (n={row[0]})"))
            handles.append(handle)
        ax.legend(title = 'Cell condition', handles=handles, bbox_to_anchor=bbox)
        if save_pdf == 'No':
            plt.savefig(os.path.join(outdir, f'{chip}_{metric}_boxplots_all_{facet_by}'), bbox_inches = 'tight')
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')        
              

def _make_boxplots_for_each_pool(dataframe, metric, outdir, coloured_by = 'cell_condition', save_pdf= 'No'):
    '''makes individual boxplots for every pool thats in the dataframe output by merge_sequencing_csvs_with_metrics_csvs()
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want boxplots for 
    outdir: where you want your plots deposited 
    coloured_by: boxplots are automatically coloured by cell condition but you can change this if you want.
    '''
    poolid_list = list(dataframe.pool_id.unique())
    for poolid in poolid_list:
        poolid_subset = dataframe[dataframe["pool_id"] == poolid]
        plt.figure()
        sns.boxplot(data=poolid_subset, x=coloured_by, y=metric).set(
            yscale = "log",
            xlabel = make_metric_name_look_nice(coloured_by), 
            ylabel = make_metric_name_look_nice(metric), 
            title = str(poolid + ' ' + make_metric_name_look_nice(metric)))
        if save_pdf == 'No':
            plt.savefig(os.path.join(outdir, str(poolid + '_' + metric + "_boxplot.png")))
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')

def make_legend_handles(colour_dict):
    '''Needs a colour dictionary with keys as the condition you want the legend for, and values as colours'''
    handle_list = []
    for condition in colour_dict.keys():
        condition_leg = mlines.Line2D([], [], color=colour_dict[condition], markersize=15, label=condition)
        handle_list.append(condition_leg)

def make_simple_boxplot(dataframe, xmetric, ymetric,  outdir, log = False, save_pdf= 'No', custom_legend = 'No', palette = 'default', hue = 'None', bbox = (1,1)):
    '''makes a simple boxplot from an x and y metric
    outdir: where to save plots
    log: set to True if you want a log y scale
    save_pdf: used to save in a pdf with other plots
    custom_legend: if you want a custom legend, feed this handles produced by make_legend_handles'''
    df_num_cells = dataframe.drop_duplicates().groupby(xmetric).size().reset_index()
    unique_x_conditions = dataframe[xmetric].unique()
    if palette == 'default':
        pal = sns.color_palette()
    else:
        pal = palette
    plt.figure()
    if hue == 'None':
        lplot = sns.boxplot(data=dataframe, x=xmetric, y=ymetric, palette = pal).set(
            xlabel = make_metric_name_look_nice(xmetric), 
            ylabel = make_metric_name_look_nice(ymetric), 
            title = f"{make_metric_name_look_nice(xmetric)} vs. {make_metric_name_look_nice(ymetric)}")
    else:
        lplot = sns.boxplot(data=dataframe, x=xmetric, y=ymetric, hue = hue, palette = pal, dodge = False).set(
            xlabel = make_metric_name_look_nice(xmetric), 
            ylabel = make_metric_name_look_nice(ymetric), 
            title = f"{make_metric_name_look_nice(xmetric)} vs. {make_metric_name_look_nice(ymetric)}")      
    plt.xticks(rotation = 45)
    if log == True:
        plt.yscale('log')
        
    if custom_legend != 'No':
        lplot.legend(handles=[custom_legend], bbox_to_anchor=bbox)
    else: 
        plt.legend(bbox_to_anchor=bbox)
    if save_pdf == 'No':
        plt.savefig(os.path.join(outdir, f"{xmetric}_{ymetric}_boxplot.png"), bbox_inches = 'tight')
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')

def make_density_lineplot(dataframe, outdir, xmetric, ymetric, colour_by = 'cell_condition', palette = colour_dict, split_by = 'cell_id', save_pdf='No', log = False, bbox = (1,1)):
    '''makes a lineplot meant to represent the distribution of a given metric for all cells in dataframe 
   dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    xmetric: the metric the ditribution is displayed over (xaxis)
    ymetric: the number of observations that fall in a particular spot in the distribution (yaxis)
    outdir: where you want the figure to be saved 
    colour_by: variable to split the lineplot colouring by - if you change this you must also change palette argument
    palette: colour palette dictionary used to assign colour to variable specified by 'colour_by'''
    dataframe = dataframe.dropna(subset=[xmetric, ymetric])
    dataframe = dataframe.astype({xmetric: 'int', ymetric:'int'})
    chipid_list = list(dataframe.chip_id.unique())
    for chip in chipid_list:
        dataframe = dataframe.loc[dataframe['chip_id'] == chip] 
        fig, ax = plt.subplots()
        plt.title(str(str(chip) + " " + make_metric_name_look_nice(xmetric) + " density lineplot"))
        plt.xticks(rotation = 45)
        if split_by == 'cell_id':
            sns.lineplot(data=dataframe, x=xmetric, y=ymetric, hue=colour_by, palette = palette, units=split_by, estimator=None, alpha=0.1).set(
                xlabel = make_metric_name_look_nice(xmetric),
                ylabel = make_metric_name_look_nice(ymetric)
            )
        else:
            sns.lineplot(data=dataframe, x=xmetric, y=ymetric, hue=colour_by, palette = palette).set(
                xlabel = make_metric_name_look_nice(xmetric),
                ylabel = make_metric_name_look_nice(ymetric)
            )
        ax.legend(title = 'Cell condition', bbox_to_anchor=bbox)
        if log == True:
            plt.yscale('log')
            
        if save_pdf == 'No':
            plt.savefig(os.path.join(outdir, str(str(chip) + '_' + xmetric + '_density_lineplot.png')), bbox_inches = 'tight')
            plt.close('all')
        else: 
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all') 


def _make_multi_axis_lineplot(dataframe, xmetric, y1metric, y2metric, outdir, plot_all_cells = "all_cells", set_y_lim = "No", save_pdf='No'):
    '''
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    Makes line plots with multiple y axis
    plot_all_cells: determines if 'units' argument in seaborn will be used to make a line for each and every cell. 
        this should be set to either the Automatic "all_cells" or to "cell_conditions" to give the plots proper titles.
    set_y_lim: some metrics (gc normalized coverage) need a y lim set. Default is no limit. Any other argument must be numberical value that will be set as the upper y axis limit
    '''
    for chip_id in dataframe.chip_id.unique():
        
        if plot_all_cells == "all_cells":
            ax1 = sns.lineplot(data=dataframe, x=xmetric, y=y1metric, hue='cell_condition', palette = colour_dict, units='cell_id', estimator=None, alpha=0.01)
        else:
            ax1 = sns.lineplot(data=dataframe, x=xmetric, y=y1metric, hue='cell_condition', palette = colour_dict, alpha=0.5)
        if set_y_lim != 'No':
            ax1.set_ylim(0,set_y_lim)
        ax1.set_ylabel(make_metric_name_look_nice(y1metric))
        ax1.set_xlabel(make_metric_name_look_nice(xmetric))
        ax1.legend(loc = "upper left")

        ax2 = ax1.twinx()
        ax2 = sns.lineplot(data=dataframe, x=xmetric, y=y2metric, color='black')
        ax2.set_ylabel(make_metric_name_look_nice(y2metric))
        ax2.legend([make_metric_name_look_nice(y2metric)], loc = "upper right")
        title = str(chip_id + f" {make_metric_name_look_nice(xmetric)} {make_metric_name_look_nice(y1metric)} for {make_metric_name_look_nice(plot_all_cells)}")
        plt.title(title, y=1.05)
        plt.xticks(rotation = 45)

        if save_pdf=='No':
            plt.savefig(os.path.join(outdir, f"{chip_id}_{xmetric}_{y1metric}_{plot_all_cells}_lineplot.png"))
            plt.close('all')
        else: 
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')

def make_multiple_lineplots_with_same_ax(dataframe, xmetric, y1metric, y2metric, outdir, plot_all_cells = "all_cells", colour_by = 'cell_condition', palette = colour_dict, set_y_lim = "No", save_pdf='No'):
    for chip_id in dataframe.chip_id.unique():
        fig, ax = plt.subplots()
        
        if palette == 'automatic':
            palette = make_colour_dict(dataframe.colour_by.unique())
        else:
            palette = palette 
            
        if plot_all_cells == "all_cells":
            y1_plot = sns.lineplot(data=dataframe, x=xmetric, y=y1metric, hue=colour_by, palette = palette, units='cell_id', estimator=None, alpha=0.01, ax=ax)
            title = str(f"{chip_id} {make_metric_name_look_nice(y1metric)} for {make_metric_name_look_nice(plot_all_cells)}")
        else:
            y2_plot = sns.lineplot(data=dataframe, x=xmetric, y=y1metric, hue=colour_by, palette = palette, alpha=0.5, ax=ax)
            title = str(f"{chip_id} {make_metric_name_look_nice(y1metric)} for Aggregated Cells")
        plt.title(title,y=1.05)
        
        if set_y_lim != 'No':
            ax.set_ylim(0,set_y_lim)
        ax.set_ylabel(make_metric_name_look_nice(y1metric))
        ax.set_xlabel(make_metric_name_look_nice(xmetric))
        y1_leg = plt.legend(loc=('upper left'), title='Cell Condition')

        plt.gca().add_artist(y1_leg)

        sns.lineplot(data=dataframe, x=xmetric, y=y2metric, color='black', ax=ax)
        black_line = mlines.Line2D([], [], color='black',
                            markersize=15, label='Reference GC Content')
        plt.legend(handles=[black_line], loc=('upper right'))
        
        
        if y1metric == 'normalized_coverage':
            ax.axhline(y=1.0, color='m')

        if save_pdf=='No':
            plt.savefig(os.path.join(outdir, f"{chip_id}_gc_plot.png"))
            plt.close('all')
        else: 
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')



def make_scatter(dataframe, outdir, xmetric, ymetric, title = 'default', colour_by = 'cell_condition', save_pdf='No'):
    '''makes a scatter plot of xmetric by ymetric an ploops it into the specified outdir
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() '''
    chipid_list = list(dataframe.chip_id.unique())
    for chipid in chipid_list:
        dataframe = dataframe.loc[dataframe['chip_id'] == chipid] 
        plt.figure()
        ax = sns.scatterplot(data=dataframe, x=xmetric, y=ymetric, hue = colour_by, palette = colour_dict).set(
            xlabel = make_metric_name_look_nice(xmetric), 
            ylabel = make_metric_name_look_nice(ymetric)
        )

        if title != 'defualt':
            plt.title(title)
        
        if save_pdf=='No':
            plt.savefig(os.path.join(outdir, str(chipid + '_' + xmetric + '_' + ymetric + '_scatterplot.png')), bbox_inches = 'tight')
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')


def make_plate_matrix(size = 72):
    '''
    Makes a dataframe the size of the isolatrix plate (72x72).
    You need this in order to generate a heatmap of the right size - you merge this with your df to make it the right size '''
    df = pd.DataFrame()
    for row in range(0,size):
        column_list = []
        for column in range(0,size):
            col_value = (1 + column,1 + row)
            column_list.append(col_value)
        column_array = pd.array(column_list)
        df[str(row + 1)] = column_array.tolist()
    df['row'] = range(1, size+1)
    df['column'] = range(1, size+1)
    return df


def make_heatmap(dataframe, metric, outdir, plate_size=72, save_pdf='No'):
    '''
    make a heatmap for each chip in a dataframe that represents the value of a given metric in each well.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want to make a heatmap out of
    outdir: where you want the heatmap to be saved 
    plate size: is the size of the isolatrix plate - 72x72
    '''
    plate_matrix = make_plate_matrix(plate_size)
    chipid_list = list(dataframe.chip_id.unique())
    for chipid in chipid_list:
        chipid_subset = dataframe[dataframe["chip_id"] == chipid]
        df_merge = pd.merge(chipid_subset, plate_matrix, how = "outer", on=['row', 'column'])
        df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
        sns.set_theme(style='whitegrid')
        sns.set_context("talk")
        plt.figure()
        fig, ax = plt.subplots(figsize=(13,7))
        title = (str(chipid) + ' ' + str(make_metric_name_look_nice(metric)))
        plt.title(title)
        ttl = ax.title
        ttl.set_position([0.5,1.05])
        sns.heatmap(df_pivot,linewidths=0.30,cmap = "flare", ax=ax)
        sns.reset_defaults()
        sns.set_context("talk")
        if save_pdf == 'No':
            plt.savefig(os.path.join(outdir, (str(chipid) + "_" + str(metric) + "_heatmap.png")))
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')
    
#make_heatmap(merged_df, 'number_of_read_pairs', out_dir)



def make_chipid_cellcondition_heatmaps(dataframe, metric, outdir, plate_size=72, save_pdf='No'):
    ''''
    Makes a different heatmap for each cell condition for a given metric on a chip.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want to make a heatmap out of
    outdir: where you want the heatmap to be saved 
    plate size: is the size of the isolatrix plate - 72x72
    '''
    plate_matrix = make_plate_matrix(plate_size)
    chipid_list = list(dataframe.chip_id.unique())
    for chipid in chipid_list:
        chipid_subset = dataframe[dataframe["chip_id"] == chipid]
        max_metric_val = chipid_subset[metric].max()
        min_metric_val = chipid_subset[metric].min()
        df_merge = pd.merge(chipid_subset, plate_matrix, how = "outer", on=['row', 'column'])
        df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
        sns.set_theme(style='whitegrid')
        sns.set_context("talk")
        plt.figure()
        fig, ax = plt.subplots(figsize=(13,7))
        title = (str(chipid) + ' ' + str(make_metric_name_look_nice(metric)))
        plt.title(title, fontsize='large')
        ttl = ax.title
        ttl.set_position([0.5,1.05])
        sns.heatmap(df_pivot,linewidths=0.30,cmap = "flare", ax=ax, vmin = min_metric_val, vmax=max_metric_val)
        if save_pdf == "No":
            plt.savefig(os.path.join(outdir, (str(chipid) + "_" + str(metric) + "_heatmap.png")))
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')
        cell_condition_list = list(chipid_subset.cell_condition.unique())
        for cell_condition in cell_condition_list:
            condition_subset = chipid_subset[chipid_subset["cell_condition"] == cell_condition]
            condition_merge = pd.merge(condition_subset, plate_matrix, how = "outer")
            condition_pivot = condition_merge.pivot(index='row',columns='column',values=metric)
            plt.figure()
            fig, ax = plt.subplots(figsize=(13,7))
            title = str(chipid + " " + make_metric_name_look_nice(cell_condition) + ' ' + make_metric_name_look_nice(metric))
            plt.title(title, fontsize='large')
            ttl = ax.title
            ttl.set_position([0.5,1.05])
            sns.heatmap(condition_pivot,linewidths=0.30,cmap = "flare", ax=ax, vmin = min_metric_val, vmax=max_metric_val)
            sns.reset_defaults()
            sns.set_context("talk")
            if save_pdf == 'No':
                plt.savefig(os.path.join(outdir, str(chipid + "_" + cell_condition + "_" + metric + "_heatmap.png")))
                plt.close('all')
            else:
                save_pdf.savefig(bbox_inches = 'tight')
                plt.close('all')

def make_chipid_cellcondition_heatmaps_no_fixed_scale(dataframe, metric, outdir, plate_size=72, save_pdf='No'):
    ''''
    Makes a different heatmap for each cell condition for a given metric on a chip.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want to make a heatmap out of
    outdir: where you want the heatmap to be saved 
    plate size: is the size of the isolatrix plate - 72x72
    '''
    plate_matrix = make_plate_matrix(plate_size)
    chipid_list = list(dataframe.chip_id.unique())
    for chipid in chipid_list:
        chipid_subset = dataframe[dataframe["chip_id"] == chipid]
        max_metric_val = chipid_subset[metric].max()
        min_metric_val = chipid_subset[metric].min()
        df_merge = pd.merge(chipid_subset, plate_matrix, how = "outer", on=['row', 'column'])
        df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
        sns.set_theme(style='whitegrid')
        sns.set_context("talk")
        plt.figure()
        fig, ax = plt.subplots(figsize=(13,7))
        title = (str(chipid) + ' ' + str(make_metric_name_look_nice(metric)) + ' variable scale')
        plt.title(title, fontsize='large')
        ttl = ax.title
        ttl.set_position([0.5,1.05])
        sns.heatmap(df_pivot,linewidths=0.30,cmap = "flare", ax=ax)
        if save_pdf == "No":
            plt.savefig(os.path.join(outdir, (str(chipid) + "_" + str(metric) + "_heatmap.png")))
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')
        cell_condition_list = list(chipid_subset.cell_condition.unique())
        for cell_condition in cell_condition_list:
            condition_subset = chipid_subset[chipid_subset["cell_condition"] == cell_condition]
            condition_merge = pd.merge(condition_subset, plate_matrix, how = "outer")
            condition_pivot = condition_merge.pivot(index='row',columns='column',values=metric)
            plt.figure()
            fig, ax = plt.subplots(figsize=(13,7))
            title = str(chipid + " " + make_metric_name_look_nice(cell_condition) + ' ' + make_metric_name_look_nice(metric) + ' variable scale')
            plt.title(title, fontsize='large')
            ttl = ax.title
            ttl.set_position([0.5,1.05])
            sns.heatmap(condition_pivot,linewidths=0.30,cmap = "flare", ax=ax)
            sns.reset_defaults()
            sns.set_context("talk")
            if save_pdf == 'No':
                plt.savefig(os.path.join(outdir, str(chipid + "_" + cell_condition + "_" + metric + "_heatmap.png")))
                plt.close('all')
            else:
                save_pdf.savefig(bbox_inches = 'tight')
                plt.close('all')



def make_contam_boxplot(fastq_screen_csv, iso_csv, save_pdf= 'No', make_all = True, bbox = (1.8,1)):
    '''makes a boxplot displaying collated output from fastqscreen to indicate contamination in run
    fastq_screen_csv: str, is the path to the csv output by collate_fastq_screen function
    iso_csv: str, path to the csv given for the sequencing run '''
    fastq_screen_df = pd.read_csv(fastq_screen_csv)
    iso_df = pd.read_csv(iso_csv)
    iso_df['cell_id'] = (iso_df['chip_id'] + '-' + iso_df['pool_id'] + '-R' + iso_df['row'].apply(str) + '-C' + iso_df['column'].apply(str))
    merged = fastq_screen_df.merge(iso_df, how = 'left', on = 'cell_id')
    merged = merged.loc[merged['variable'].isin(['one_hit_one_genome', 'multiple_hits_one_genome'])]
    merged['number_of_reads'] = merged['value']

    chipid_list = list(merged.chip_id.unique())
    for chip_id in chipid_list:
        fig, ax = plt.subplots()
        sns.barplot(data = merged, x = 'genome', y = 'number_of_reads', hue = 'variable').set(
                                                                                              ylabel = 'Number of reads ',
                                                                                              xlabel = 'Genome')
        plt.xticks(rotation = 90)
        plt.title((str(chip_id) + ' Number of reads mapping to different genomes'), y=1.1)
        ax.legend().set_title('')
        one_hit_leg = mlines.Line2D([], [], color='#1f77b4',
                                markersize=15, label='One hit one genome')
        multi_hit_leg = mlines.Line2D([], [], color='#ff7f0e',
                        markersize=15, label='Multiple hits one genome')
        ax.legend(handles=[one_hit_leg, multi_hit_leg], bbox_to_anchor=bbox)
        
        print(chip_id)
        if save_pdf == 'No':
            plt.savefig(os.path.join(outdir, str(chip_id + "_contam_boxplot.png")))
            plt.close('all')
        else:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')

        if make_all == True: 
            for cell_condition in merged.cell_condition.unique():
                merged_subset = merged[merged['cell_condition'] == cell_condition]
                fig, ax = plt.subplots()
                sns.barplot(data = merged_subset, x = 'genome', y = 'number_of_reads', hue = 'variable').set(
                                                                                                    ylabel = 'Number of reads ',
                                                                                                    xlabel = 'Genome')
                #plt.yscale('log')
                plt.title((str(chip_id) + ' ' + str(cell_condition) + ' Number of reads mapping to different genomes'), y=1.1)
                plt.xticks(rotation = 90)
                ax.set_ylim([0, 3000]) 
                ax.legend().set_title('')
                one_hit_leg = mlines.Line2D([], [], color='#1f77b4',
                                        markersize=15, label='One hit one genome')
                multi_hit_leg = mlines.Line2D([], [], color='#ff7f0e',
                                markersize=15, label='Multiple hits one genome')
                ax.legend(handles=[one_hit_leg, multi_hit_leg], loc=('upper right'))
                
                
                if save_pdf == 'No':
                    plt.savefig(os.path.join(outdir, str(chip_id + "_contam_boxplot.png")))
                    plt.close('all')
                else:
                    save_pdf.savefig(bbox_inches = 'tight')
                    plt.close('all')


def make_breadth_depth_boxplot(df, outdir, save_pdf='No', bbox = (1,1)):
    '''df = collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip 
outdir = where you want figure to be saved
bbox = adjusts legend position  
'''
    chip = df.chip_id.unique()[0]
    df = df[['cell_id', 'cell_condition', 'mean_coverage', 'coverage_breadth']]
    df = df.rename({'mean_coverage': 'mean_coverage_depth'}, axis='columns')
    df = df.melt(id_vars=['cell_id', 'cell_condition'], value_vars=['mean_coverage_depth', 'coverage_breadth'])
    fig,ax = plt.subplots()
    sns.boxplot(x='cell_condition', y='value', hue = 'variable', data=df).set(
        ylabel = ' ',
        xlabel = 'Cell condition',
        title = f'{chip} Mean coverage depth and breadth for each cell condition'
    )
    mean_cov = mlines.Line2D([], [], color='#1f77b4',
                            markersize=15, label='Mean coverage depth')
    depth = mlines.Line2D([], [], color='#ff7f0e',
                    markersize=15, label='Coverage breadth')
    ax.legend(handles=[mean_cov, depth], bbox_to_anchor=bbox)     
    if save_pdf == 'No':   
        plt.savefig(os.path.join(outdir, 'breadth_depth_boxplot.png'), bbox_inches = 'tight')
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')
    


def make_total_raw_reads_plots(dataframe, outdir, save_pdf='No'):
        _make_faceted_combined_boxplots(dataframe, 'total_raw_reads', outdir, save_pdf=save_pdf, log='yes')
        _make_boxplots_for_each_pool(dataframe, 'total_raw_reads', outdir, save_pdf=save_pdf)
        make_heatmap(dataframe, 'total_raw_reads', outdir, save_pdf=save_pdf)
        _make_chipid_cellcondition_heatmaps(dataframe, 'total_raw_reads', outdir, save_pdf=save_pdf)


def make_insert_hist_df(collated_input_dir, df):
    ''' Makes a dataframe for the insert size histogram that can be used for insert distribution plotting 
collated_input_dir: directory where you have your collated output from the bam_metrics file
df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip '''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id', 'pool_id']]
    chip_id = cell_condition_df.chip_id[0]
    for file in os.listdir(collated_input_dir):
        if str(chip_id + '.insert_hist') in file: 
            df_insert = pd.read_csv(os.path.join(collated_input_dir, file))
            df_insert = pd.melt(df_insert, id_vars=['cell_id'], var_name='insert_size', value_name='number_of_reads')
            df_insert.cell_id.astype(str)
            df_insert['insert_size'] = df_insert.insert_size.astype(float)
            df_insert = df_insert.merge(cell_condition_df, how='inner', on='cell_id')
            #exclude outliers with massive insert sizes
            df_insert = df_insert[df_insert['insert_size'] <= 2000]
            return df_insert

        

def make_insert_cell_condition_faceted_lineplot(dataframe, outdir, facet_by = 'cell_condition', colour_dict = 'automatic', save_pdf = 'No'):
    ''' Makes lineplots for each cell condition in one large faceted plot 
    df: insert histogram df for a chip produced by _make_insert_hist_df
    outdir: where to save resulting plots 
    facet_by: the metric you want to facet the plots by 
    colour_dict: if set to automatic will automaticlaly generate a colour dict based on seaborn default colours, other option is to set it to a predefined dictionary 
    '''
    facet_condition_list = list(dataframe[facet_by].unique())
    if colour_dict == 'automatic':
        colour_dict = make_colour_dict(facet_condition_list)
    else:
        colour_dict = colour_dict
    chip_id = dataframe.chip_id[0]
    fig, axes = plt.subplots(nrows=len(facet_condition_list), ncols=1, figsize=(8, 3*len(facet_condition_list)), sharex=True, sharey=False)
    if len(facet_condition_list) == 4:
        suptitle_y = 0.96
    elif len(facet_condition_list) == 3:
        suptitle_y = 0.98
    elif len(facet_condition_list) == 2:
        suptitle_y = 1.05
    elif len(facet_condition_list) == 1:
        suptitle_y = 1.2
    else: 
        suptitle_y = 0.96
    
    fig.suptitle('Insert size distribution for all cells', y = suptitle_y, size = 18)
    #fig.suptitle(f"""{chip_id} {make_metric_name_look_nice(facet_by)} 
#insert size distribution for all cells""", y = suptitle_y, size = 18)
    max_reads = dataframe.number_of_reads.max()
    for i, condition in enumerate(facet_condition_list):
        dataframe_subset = dataframe.loc[dataframe[facet_by] == condition]
        if len(facet_condition_list) == 1:
            lplot = sns.lineplot(data=dataframe_subset, x="insert_size", y="number_of_reads", color = colour_dict[condition], units='cell_id', estimator=None, alpha=0.03)
        else: 
            lplot = sns.lineplot(ax=axes[i],data=dataframe_subset, x="insert_size", y="number_of_reads", color = colour_dict[condition], units='cell_id', estimator=None, alpha=0.03)
        lplot.set_ylabel('Number of paired reads', fontsize=14)
        lplot.set_xlabel('Insert size', fontsize = 14)
        lplot.set_ylim(0,max_reads)
        lplot.set_xlim(0,1200)
        condition_leg = mlines.Line2D([], [], color=colour_dict[condition],
                        markersize=15, label=condition)
        lplot.legend(handles=[condition_leg], loc=('upper left'))    
        fig = lplot.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
    if save_pdf == 'No':
        plt.savefig(os.path.join(outdir, "insert_size_lineplot.png"), bbox_inches = 'tight')
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')

#pivot histos to make density plots
def make_insert_plots(collated_input_dir, df, outdir, save_pdf='No', make_all = True):
    '''makes a cute little density plot of individual cell's insert size distribution coloured by their cell condition.
    collated_input_dir: directory where you have your collated output from the bam_metrics file
    df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip 
    outdir: where you want your nice plot deposited
    save_pdf: whether to save all plots in one pdf - default is to save as multiple pngs
    make_all: make all the insert plots (cells aggregated and split by condition in same plot, cells not aggregated and split by condition, )
    '''
    df_insert = make_insert_hist_df(collated_input_dir, df)
    if make_all == True:
        #make one of all cell individually plot coloured by cell condition 
        make_density_lineplot(df_insert, outdir, "insert_size", "number_of_reads", save_pdf = save_pdf)
    #make one plot of cells coloured and aggregated by cell condition 
    make_density_lineplot(df_insert, outdir, "insert_size", "number_of_reads", split_by = 'dont split', save_pdf=save_pdf, bbox = (1.4, 1))
    #make insert distribution lineplots of all cells, faceted and coloured by cell condition 
    make_insert_cell_condition_faceted_lineplot(df_insert, outdir, save_pdf = save_pdf)
      

def make_gc_hist_df(collated_input_dir, df):
    ''' Makes a dataframe for the gc content histogram that can be used for gc distribution plotting 
collated_input_dir: directory where you have your collated output from the bam_metrics file
df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip'''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id', 'pool_id']]
    chip_id = cell_condition_df.chip_id[0]
    for file in os.listdir(collated_input_dir):
        if str(chip_id +'.gc_bias_metrics') in file: 
            df_metrics = pd.read_csv(os.path.join(collated_input_dir, file))
            df_metrics = df_metrics.merge(cell_condition_df, how='inner', on='cell_id')
            df_metrics['reference_gc_content'] = df_metrics['windows'] / 100000000
            df_metrics['gc_content'] = df_metrics['gc']
    return df_metrics

def make_wgs_hist_df(collated_input_dir, df):
    ''' Makes a dataframe for the wgs histogram that can be used for uniformity plotting 
collated_input_dir: directory where you have your collated output from the bam_metrics file
df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip'''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id', 'pool_id']]
    chip_id = cell_condition_df.chip_id[0]
    for file in os.listdir(collated_input_dir):
        if str(chip_id +'.wgs_hist.csv') in file: 
            df_metrics = pd.read_csv(os.path.join(collated_input_dir, file))
            df_metrics = df_metrics.merge(cell_condition_df, how='inner', on='cell_id')
            #df_metrics['gc_content'] = df_metrics['gc']
        box_range = []
    for i in range(0,251):
        box_range.append(str(i))
    df_metrics = pd.melt(df_metrics, id_vars=['cell_id', 'cell_condition', 'chip_id', 'pool_id'], value_vars=box_range, var_name='coverage_depth', value_name='n_bases', ignore_index=False)
    df_metrics['coverage_depth'] = df_metrics['coverage_depth'].astype(int)
    print(df_metrics.coverage_depth.dtype)
    return df_metrics

def gc_dist_plot_aggregated(df_metrics, outdir, save_pdf = 'No'):
    '''Normalized coverage GC one plot - plots cells aggregated and coloured by cell condition
    df_metrics: produced by make_gc_hist_df() for a particular chip'''
    make_multiple_lineplots_with_same_ax(df_metrics, 'gc_content', 'normalized_coverage', 'reference_gc_content', outdir, set_y_lim = 2, plot_all_cells = "cell_condition", save_pdf=save_pdf)
     
     
def gc_dist_plot_single_cells(df_metrics, outdir, save_pdf = 'No', colour_by = 'cell_condition', palette = colour_dict):
    '''Normalized coverage GC one plot - plots all cells, colours by cell condition
    df_metrics: produced by make_gc_hist_df() for a particular chip'''
    make_multiple_lineplots_with_same_ax(df_metrics, 'gc_content', 'normalized_coverage', 'reference_gc_content', outdir, colour_by = colour_by, palette = palette, set_y_lim = 2, save_pdf=save_pdf)


def gc_cell_condition_faceted_lineplot(df_metrics, outdir, facet_by = 'cell_condition', colour_dict = 'automatic', save_pdf = 'No'):
    '''Normalized coverage GC faceted plot - plots all cells, coloured and faceted by cell condition
    df_metrics: produced by make_gc_hist_df() for a particular chip
    facet_by: the metric you want to facet the plots by 
    colour_dict: if set to automatic will automaticlaly generate a colour dict based on seaborn default colours, other option is to set it to a predefined dictionary '''
    dataframe = df_metrics
    chip_id = dataframe.chip_id[0]
    facet_condition_list = list(dataframe[facet_by].unique())
    fig, axes = plt.subplots(nrows=len(facet_condition_list), ncols=1, figsize=(8, 3*len(facet_condition_list)), sharex=True, sharey=False)
    
    if len(facet_condition_list) == 4:
        suptitle_y = 0.96
    elif len(facet_condition_list) == 3:
        suptitle_y = 0.98
    elif len(facet_condition_list) == 2:
        suptitle_y = 1.05
    elif len(facet_condition_list) == 1:
        suptitle_y = 1.2
    else: 
        suptitle_y = 0.96
        
    fig.suptitle(f"""{chip_id} {make_metric_name_look_nice(facet_by)} 
normalized coverage for all cells""", y = suptitle_y, size = 18)
    
    for i, condition in enumerate(facet_condition_list):
        dataframe_subset = dataframe.loc[dataframe[facet_by] == condition]
        if len(facet_condition_list) == 1:
            lplot = sns.lineplot(data=dataframe_subset, x='gc_content', y='normalized_coverage', color = colour_dict[condition], units='cell_id', estimator=None, alpha=0.03)
            refplot = sns.lineplot(data=dataframe_subset, x='gc_content', y='reference_gc_content', color = 'black')
        else: 
            lplot = sns.lineplot(ax=axes[i],data=dataframe_subset, x='gc_content', y='normalized_coverage', color = colour_dict[condition], units='cell_id', estimator=None, alpha=0.03)
            refplot = sns.lineplot(ax=axes[i], data=dataframe_subset, x='gc_content', y='reference_gc_content', color = 'black')
        
        lplot.set_ylim(0,2)
        lplot.set_ylabel('Normalized coverage', fontsize = 14)
        lplot.set_xlabel('GC content', fontsize = 14)
        lplot.axhline(y=1.0, color='m')
        condition_leg = mlines.Line2D([], [], color=colour_dict[condition],
                        markersize=15, label=condition)
        reference_leg = mlines.Line2D([], [], color='black',
                        markersize=15, label='Reference GC Content')
        lplot.legend(handles=[condition_leg, reference_leg], loc=('upper left'))
        fig = lplot.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
    
    if save_pdf == 'No':
        plt.savefig(os.path.join(outdir, str(chip_id + '_' + facet_by + '_faceted_gc_lineplot.png')), bbox_inches = 'tight')
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')
    

def make_gc_plots(collated_input_dir, df, outdir, save_pdf='No',make_all = True):
    '''makes relevant plots to visually display GC bias results from a sequencing run - meant to be used after the bam_metrics pipeline
    collated_input_dir: directory where you have your collated output from the bam_metrics file
    df: dataframe you should get from running merge_sequencing_csvs_with_metrics_csvs()
    outdir: where you want your nice plot deposited'''
    sns.reset_defaults()
    
    df_metrics = make_gc_hist_df(collated_input_dir, df)
    
    #Normalized coverage GC one plot - plots all cells, colours by cell condition
    if make_all == True:
       gc_dist_plot_single_cells(df_metrics, outdir, save_pdf = save_pdf)
    #Normalized Coverage Pooled Cell Conditions 
    gc_dist_plot_aggregated(df_metrics, outdir, save_pdf = save_pdf)

    gc_cell_condition_faceted_lineplot(df_metrics, outdir, save_pdf = save_pdf)
    
    sns.set_context("talk")      



def make_duplicate_plots(df, outdir, save_pdf = 'No'):
    '''makes relevant plots to visually display duplicate metrics from a sequencing run - meant to be used after the bam_metrics pipeline
    df: dataframe you should get from running merge_sequencing_csvs_with_metrics_csvs() for a particular chip 
    outdir: where you want your nice plot deposited'''
    _make_faceted_combined_boxplots(df, 'duplicate_fraction', outdir, log = 'no', save_pdf=save_pdf, bbox = (1.6, 1))


def make_wgs_plots(collated_input_dir, df, outdir, save_pdf = 'No', make_all = True):
    '''makes relevant plots to visually display wgs metrics results from a sequencing run - meant to be used after the bam_metrics pipeline
    collated_input_dir: directory where you have your collated output from the bam_metrics file
    df: dataframe you should get from running merge_sequencing_csvs_with_metrics_csvs()
    outdir: where you want your nice plot deposited'''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id', 'pool_id', 'row', 'column']]
    chip_id = df.chip_id.unique()[0]
    for file in os.listdir(collated_input_dir):
        if make_all == True:
            if str(chip_id + '.wgs_hist') in file: 
                df_wgs_hist = pd.read_csv(os.path.join(collated_input_dir, file))
                df_wgs_hist = pd.melt(df_wgs_hist, id_vars=['cell_id'], var_name='coverage', value_name='high_quality_coverage_count')
                df_wgs_hist = df_wgs_hist.merge(cell_condition_df, how='left', on='cell_id')
                df_wgs_hist = df_wgs_hist[df_wgs_hist['coverage'].astype('int') > 0]
                make_density_lineplot(df_wgs_hist, outdir, "coverage", "high_quality_coverage_count", save_pdf=save_pdf, log = True)
        if str(chip_id + '.wgs_metrics') in file: 
            df_wgs_metrics = pd.read_csv(os.path.join(collated_input_dir, file))
            df_wgs_metrics = df_wgs_metrics.merge(cell_condition_df, how='left', on='cell_id')
            df_wgs_metrics['mean_coverage_depth'] = df_wgs_metrics['mean_coverage']
            df_wgs_metrics.reset_index(drop=True)
            _make_chipid_cellcondition_heatmaps(df_wgs_metrics, 'mean_coverage_depth', outdir, save_pdf=save_pdf)
            if make_all == True:
                _make_chipid_cellcondition_heatmaps_no_fixed_scale(df_wgs_metrics, 'mean_coverage_depth', outdir, save_pdf=save_pdf)
                _make_faceted_combined_boxplots(df_wgs_metrics, 'mean_coverage_depth', outdir, save_pdf=save_pdf, log = False)


def make_mapping_plots(df, outdir, save_pdf, make_all = True):
    if make_all == True: 
        _make_faceted_combined_boxplots(df, 'mapping_percent', outdir, save_pdf=save_pdf, log = False)
    _make_chipid_cellcondition_heatmaps(df, 'mapping_percent', outdir, save_pdf=save_pdf)

    


def make_all_plots_for_library_quality(collated_input_dir, iso_csv, chip_id, fastq_screen_csv, outdir):
    '''
    collated_input_dir: path to directory containing collated output from bam_metrics pipeline
    iso_csv: sequencing csv that matches the isolatrix run 
    outdir: path to where you want your pdf to be saved 
    '''
    df = merge_sequencing_csvs_with_metrics_csvs(collated_input_dir, iso_csv, chip_id)
    seq_run = iso_csv.split('/')[-1].split('.')[0]
    outfile_dict = { 'read_pair' : os.path.join(outdir, f'{seq_run}_readpair_plots.pdf'),
                     'insert' : os.path.join(outdir, f'{seq_run}_insert_plots.pdf'),
                     'gc': os.path.join(outdir, f'{seq_run}_gc_plots.pdf'),
                     'duplicate': os.path.join(outdir, f'{seq_run}_duplicate_plots.pdf'),
                     'wgs': os.path.join(outdir, f'{seq_run}_wgs_plots.pdf'),
                     'mapping': os.path.join(outdir, f'{seq_run}_mapping_plots.pdf'),
                     'allplots': os.path.join(outdir, f'{seq_run}_allplots.pdf'),
                     'contam' : os.path.join(outdir, f'{seq_run}_contam_plots.pdf')}



    with PdfPages(outfile_dict['mapping']) as out_pdf:
        make_mapping_plots(df, outdir, save_pdf=out_pdf)
        
    with PdfPages(outfile_dict['read_pair']) as out_pdf:
        make_total_raw_reads_plots(df, outdir, save_pdf=out_pdf)
    
    with PdfPages(outfile_dict['insert']) as out_pdf:
        make_insert_plots(collated_input_dir, df, outdir, save_pdf=out_pdf)
    
    with PdfPages(outfile_dict['gc']) as out_pdf:
        make_gc_plots(collated_input_dir, df, outdir, save_pdf=out_pdf)

    with PdfPages(outfile_dict['duplicate']) as out_pdf:
        make_duplicate_plots(df, outdir, save_pdf = out_pdf)
    
    with PdfPages(outfile_dict['wgs']) as out_pdf:
        make_wgs_plots(collated_input_dir, df, outdir, save_pdf = out_pdf)
    
    with PdfPages(outfile_dict['allplots']) as out_pdf:
        make_total_raw_reads_plots(df, outdir, save_pdf=out_pdf)
        make_insert_plots(collated_input_dir, df, outdir, save_pdf=out_pdf)
        make_gc_plots(collated_input_dir, df, outdir, save_pdf=out_pdf)
        make_duplicate_plots(df, outdir, save_pdf = out_pdf)
        make_wgs_plots(collated_input_dir, df, outdir, save_pdf = out_pdf)
    
    #with PdfPages(outfile_dict['contam']) as out_pdf:
    #    make_contam_boxplot(fastq_screen_csv, iso_csv, save_pdf= out_pdf)
    print('Successfully made all library metrics plots!') 


def make_key_library_plots(collated_input_dir, iso_csv, chip_id, fastq_screen_csv, outdir):
    '''makes one pdf with key library metrics plots 
    collated_input_dir: directory with collated gatk metrics csvs produced by bam_metrics pipeline
    iso_csv: sequencing csv that came with isolatrix run 
    fastq_screen_csv: path to the collated fastq_screen csv produced by the fastq_metrics pipeline
    outdir: path to directory for plots to be created in'''
    df = merge_sequencing_csvs_with_metrics_csvs(collated_input_dir, iso_csv, chip_id)
    chip = df.chip_id.unique()[0]
    outfile = os.path.join(outdir, f'{chip}_QC_plots.pdf')
    with PdfPages(outfile) as out_pdf:
        make_run_metrics_table(df, outdir, save_pdf = out_pdf)
        make_breadth_depth_boxplot(df, outdir, save_pdf=out_pdf)
        make_wgs_plots(collated_input_dir, df, outdir, save_pdf = out_pdf, make_all = False)
        make_mapping_plots(df, outdir, save_pdf = out_pdf, make_all = True)
        make_insert_plots(collated_input_dir, df, outdir, save_pdf=out_pdf, make_all = False)
        make_duplicate_plots(df, outdir, save_pdf = out_pdf)
        make_gc_plots(collated_input_dir, df, outdir, save_pdf = out_pdf, make_all = False)
        make_contam_boxplot(fastq_screen_csv, iso_csv, save_pdf= out_pdf, make_all = False)
print('Successfully made key library metrics plots!') 
    
    


def main():
    args = parse_args()
    #make_all_plots_for_library_quality(args.collated_input_dir, args.iso_csv, args.collated_fastq_screen_csv, args.outdir)
    make_key_library_plots(args.collated_input_dir, args.iso_csv, args.collated_fastq_screen_csv, args.outdir)

if __name__ == '__main__':
    main()