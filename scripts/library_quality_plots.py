from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
sys.path.insert(0, '../../scripts')
from generic_plots import *
import numpy as np


def make_library_metric_df(collated_input_dir, sequencing_csv, chip_id, iso_library = True, paired_end = True):
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
                    'insert_metrics':['cell_id','median_insert_size','min_insert_size','max_insert_size','mean_insert_size','standard_deviation'],
                    'wgs_metrics':['cell_id','coverage_breadth','mean_coverage','sd_coverage','median_coverage'],
                    'mt_coverage': ['cell_id', 'mt_coverage_depth'],
                    'RNA_Metrics': ['cell_id' ,'pf_bases','pf_aligned_bases', 'ignored_reads','correct_strand_reads','incorrect_strand_reads','num_r1_transcript_strand_reads','num_r2_transcript_strand_reads','num_unexplained_reads','pct_r1_transcript_strand_reads','pct_r2_transcript_strand_reads','pct_ribosomal_bases','pct_coding_bases','pct_utr_bases','pct_intronic_bases','pct_intergenic_bases','pct_mrna_bases','pct_usable_bases','pct_correct_strand_reads','median_cv_coverage','median_5prime_bias','median_3prime_bias','median_5prime_to_3prime_bias','sample','library','read_group', 'coding_bases', 'utr_bases', 'intronic_bases']}
    iso_df = pd.read_csv(sequencing_csv) 
    
    if iso_library == True:
        #create cell_id column in isolatrix/sequencing df - can't do this for other libraries 
        iso_df['cell_id'] = (iso_df['sample_id'] + '-' + iso_df['chip_id'] + '-R' + iso_df['row'].apply(str) + '-C' + iso_df['column'].apply(str))

    iso_df_merged = iso_df
    present_metrics_files = []
    i = 0
    for file in os.listdir(collated_input_dir):
        chip = file.split('.')[0]
        metric_class = file.split('.')[1]
        if chip == chip_id:
            if metric_class in desired_metrics.keys():
                try:
                    present_metrics_files.append(metric_class)
                    df = pd.read_csv(os.path.join(collated_input_dir, file))
                    df = df.loc[:,desired_metrics[metric_class]]
                    iso_df_merged = iso_df_merged.merge(df, on='cell_id', how='left')
                except:     
                    print(f"{metric_class} file is empty")

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
    iso_df_merged.replace('OddCell', 'LiveCell', inplace = True)
    #add print chip to the chip dict library for this sequencing run 
    return iso_df_merged





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


def _heatmap_code(df, metric, chipid, vmin, vmax, plate_size = 72,  title = 'default'):
    sns.set_context("talk")
    fig, ax = plt.subplots(figsize=(13,10))
    if title == 'default':
        title = f"{chipid} {str(make_metric_name_look_nice(metric))}"
    plt.title(title, fontsize='large')
    ttl = ax.title
    ttl.set_position([0.5,1.05])
    heatmap = sns.heatmap(df,linewidths=0.30, cmap = "flare", ax=ax, vmin = vmin, vmax=vmax)
    # Customize x-axis and y-axis ticks to show markers every 5 points
    xticks = np.arange(0, plate_size, 2)
    yticks = np.arange(0, plate_size, 2)
    heatmap.set_xticks(xticks)
    heatmap.set_yticks(yticks)
    # Set tick labels for both x and y axes
    heatmap.set_xticklabels(xticks, rotation=0, size = 'x-small')
    heatmap.set_yticklabels(yticks, rotation=0, size = 'x-small')


def make_heatmap(df, metric, eval_scale = True, plate_size=72, save_pdf=False, save_png = 'outplot.png', title = 'default'):
    '''
    make a heatmap for each chip in a dataframe that represents the value of a given metric in each well.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want to make a heatmap out of
    outdir: where you want the heatmap to be saved 
    plate size: is the size of the isolatrix plate - 72x72
    '''
    plate_matrix = make_plate_matrix(plate_size)
    chipid = df.chip_id.unique()[0]
    #if eval_scale == True:
    vmax = df[metric].max()
    vmin = df[metric].min()
    df_merge = pd.merge(df, plate_matrix, how = "outer", on=['row', 'column'])
    df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
    
    _heatmap_code(df_pivot, metric, chipid, vmin, vmax, plate_size = 72, title = title)

    if save_png != False:
        plt.savefig(save_png, dpi = 100, bbox_inches = 'tight')
        plt.close('all')
    if save_pdf != False:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')


def _heatmap_subplot(df, metric, vmin, vmax, ax = 0, title = ' ', change_scale = False, plate_size = 72):
    '''used to make a heatmap in one panel of a plt.subplot(), don't call until after plt.subplto has been called
        ax: the axis you want to plot into 
        axes: the axes object created by calling fig,axes = plt.subplots() - this needs to be called before running this functio
        chnage_scale = true: chnage scale to reflect each indiviudal condition you're plotting, false: keep the scale consistent to how it it for the whole chip'''
    plate_matrix = make_plate_matrix(plate_size)
    if change_scale == True:
        vmax = df[metric].max()
        vmin = df[metric].min()
    df_merge = pd.merge(df, plate_matrix, how = "outer", on=['row', 'column'])
    df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
    sns.heatmap(df_pivot,linewidths=0.30,cmap = "flare", vmin = vmin, vmax=vmax, ax=ax)
    ax.set_title(title)
    ax.set_ylabel('row')
    ax.set_xlabel('column')


def make_cellcondition_subplot_heatmaps(df, metric, outdir, plate_size=72, save_pdf='No', change_scale = False, title = ' '):
    ''' Makes a different heatmap for each cell condition for a given metric on a chip.
        dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
        metric: the metric you want to make a heatmap out of
        outdir: where you want the heatmap to be saved 
        plate size: is the size of the isolatrix plate - 72x72'''
    vmax = df[metric].max()
    vmin = df[metric].min()
    cell_conditions = df.cell_condition.unique()
    chip_id = df.chip_id.unique()[0]
    sns.set_context("paper")
    fig, axes = plt.subplots(2,2, figsize=(13,10))
    plt.suptitle(title)
    i = 0
    for condition in cell_conditions:
        condition_subset = df[df["cell_condition"] == condition]
        title = f"{condition}"
        _heatmap_subplot(condition_subset, metric, vmin, vmax, axes, ax=axes[i // 2, i % 2], title = title, change_scale = change_scale, plate_size = 72)
        i += 1
    plt.subplots_adjust(wspace=0.2, hspace=0.2, top=0.90)
    sns.reset_defaults()

    if save_pdf == 'No':
        plt.savefig(os.path.join(outdir, f"{chip_id}_ metric_heatmap_cell_condition.png"))
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')
        

def make_chipid_cellcondition_heatmaps_individually(df, metric, plate_size=72, save_pdf=False,  save_png = True, scale = 'fixed'):
    ''''
    Makes a different heatmap for each cell condition for a given metric on a chip.
    dataframe: should be a df from the dictionary produced by merge_sequencing_csvs_with_metrics_csvs() 
    metric: the metric you want to make a heatmap out of
    outdir: where you want the heatmap to be saved 
    plate size: is the size of the isolatrix plate - 72x72
    '''
    sns.set_context("paper")
    vmax = df[metric].max()
    vmin = df[metric].min()
    chipid = df.chip_id.unique()[0]

    cell_conditions = df.cell_condition.unique()
    for condition in cell_conditions:
        condition_subset = df[df["cell_condition"] == condition]
        title = f"{chipid} {condition} {make_metric_name_look_nice(metric)}"
        if scale != 'fixed':
            vmax = condition_subset[metric].max()
            vmin = condition_subset[metric].min()
            
        plate_matrix = make_plate_matrix(plate_size)
        df_merge = pd.merge(condition_subset, plate_matrix, how = "outer", on=['row', 'column'])
        df_pivot = df_merge.pivot(index='row',columns='column',values=metric)
        _heatmap_code(df_pivot, metric, chipid, vmin, vmax, plate_size = plate_size,  title = title)
        
        if save_png != False:
            plt.savefig(f"{chipid}_{condition}_{metric}_heatmap.png", bbox_inches = 'tight', dpi = 600, )
            plt.close('all')
        if save_pdf != False:
            save_pdf.savefig(bbox_inches = 'tight')
            plt.close('all')
                

def make_insert_hist_df(collated_input_dir, df):
    ''' Makes a dataframe for the experimental_conditioninsert size histogram that can be used for insert distribution plotting 
collated_input_dir: directory where you have your collated output from the bam_metrics file
df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip '''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id']]
    chip_id = f"{cell_condition_df.chip_id.unique()[0]}"
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
        


def make_insert_cell_condition_faceted_lineplot(collated_input_dir, df, facet_by = 'cell_condition', palette = 'automatic', save_pdf = False, save_png = False, title = 'default', chip_id = ' '):
    ''' Makes lineplots for each cell condition in one large faceted plot 
    df: insert histogram df for a chip produced by _make_insert_hist_df
    outdir: where to save resulting plots 
    facet_by: the metric you want to facet the plots by 
    colour_dict: if set to automatic will automaticlaly generate a colour dict based on seaborn default colours, other option is to set it to a predefined dictionary 
    '''
    dataframe = make_insert_hist_df(collated_input_dir, df)
    facet_condition_list = list(dataframe[facet_by].unique())
    if palette == 'automatic':
        pal = make_colour_dict(facet_condition_list)
    else:
        pal = palette
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
    if title == 'default':
        title = f"{chip_id} insert size distribution"
    else:
        title = title
    fig.suptitle(title, y = suptitle_y, size = 'x-large')
    max_reads = dataframe.number_of_reads.max()
    for i, condition in enumerate(facet_condition_list):
        dataframe_subset = dataframe.loc[dataframe[facet_by] == condition]
        if len(facet_condition_list) == 1:
            lplot = sns.lineplot(data=dataframe_subset, x="insert_size", y="number_of_reads", color = pal[condition], units='cell_id', estimator=None, alpha=0.03)
        else: 
            lplot = sns.lineplot(ax=axes[i],data=dataframe_subset, x="insert_size", y="number_of_reads", color = pal[condition], units='cell_id', estimator=None, alpha=0.03)
        lplot.set_ylabel('Number of paired reads', fontsize=14)
        lplot.set_xlabel('Insert size', fontsize = 14)
        lplot.set_ylim(0,max_reads)
        lplot.set_xlim(0,1200)
        condition_leg = mlines.Line2D([], [], color=pal[condition],
                        markersize=15, label=condition)
        lplot.legend(handles=[condition_leg], loc=('upper left'))    
        fig = lplot.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
    if save_png != False:
        plt.savefig(save_png, bbox_inches = 'tight', dpi = 100, )
        plt.close('all')
    if save_pdf != False:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')
        
        
def make_gc_hist_df(collated_input_dir, df):
    ''' Makes a dataframe for the gc content histogram that can be used for gc distribution plotting 
collated_input_dir: directory where you have your collated output from the bam_metrics file
df: collated_metric_df produced by merge_sequencing_csvs_with_metrics_csvs() for a particular chip'''
    cell_condition_df = df.loc[:,['cell_id', 'cell_condition', 'chip_id']]
    chip_id = f"{cell_condition_df.chip_id.unique()[0]}"
    for file in os.listdir(collated_input_dir):
        if str(chip_id +'.gc_bias_metrics') in file: 
            df_metrics = pd.read_csv(os.path.join(collated_input_dir, file))
            df_metrics = df_metrics.merge(cell_condition_df, how='inner', on='cell_id')
            df_metrics['reference_gc_content'] = df_metrics['windows'] / 100000000
            df_metrics['gc_content'] = df_metrics['gc']
    return df_metrics


def one_panel_gc_plot(gc_df, palette='automatic', aggregate_cells=True, save_png=False, save_pdf=False):
    chip_id = gc_df.chip_id.unique()[0]
    if palette == 'automatic':
        palette = make_colour_dict(gc_df.colour_by.unique())
    else:
        palette = palette
    fig, ax = plt.subplots()

    if aggregate_cells == False:
        y1_plot = sns.lineplot(data=gc_df, x='gc_content', y='normalized_coverage', hue="cell_condition", palette=palette, units='cell_id', estimator=None, alpha=0.01, ax=ax)
        title = str(f"{chip_id} normalized gc content coverage for all cells")
    else:
        y2_plot = sns.lineplot(data=gc_df, x='gc_content', y='normalized_coverage', hue="cell_condition", palette=palette, alpha=0.5, ax=ax)
        title = str(f"{chip_id} GC content coverage")

    sns.lineplot(data=gc_df, x='gc_content', y='reference_gc_content', color='black', ax=ax)
    ax.axhline(y=1.0, color='m')
    ax.set_ylim(0,2)

    plt.title(title, y=1.05)
    ax.set_ylabel('Normalized coverage')
    ax.set_xlabel('GC content')
    y1_leg = plt.legend(title='Cell Condition', bbox_to_anchor=(1, 1))

    # Create a legend for the black line
    black_line = mlines.Line2D([], [], color='black', markersize=15, label='Reference GC Content')
    
    # Get the existing handles and labels from the legend
    handles, labels = ax.get_legend_handles_labels()
    
    # Extend the handles and labels with the black line
    handles.append(black_line)
    labels.append('Reference GC Content')
    
    # Create a new legend with the combined handles and labels
    combined_leg = plt.legend(handles=handles, labels=labels, bbox_to_anchor = (1,1))
    
    # Add the combined legend to the plot
    plt.gca().add_artist(combined_leg)

    if save_png != False:
        plt.savefig(save_png, dpi=600, bbox_inches='tight')
        plt.close('all')
    if save_pdf != False:
        save_pdf.savefig(bbox_inches='tight')
        plt.close('all')




def gc_cell_condition_faceted_lineplot(df_metrics, facet_by='cell_condition', palette='automatic', save_pdf=False, save_png=False, title='default'):
    chip_id = df_metrics['chip_id'].unique()[0]
    facet_condition_list = df_metrics[facet_by].unique()
    if palette == 'automatic':
        pal = make_colour_dict(facet_condition_list)
    else:
        pal = palette
    fig, axes = plt.subplots(nrows=len(facet_condition_list), ncols=1, figsize=(8, 3 * len(facet_condition_list)), sharex=True, sharey=False)
    if len(facet_condition_list) == 4:
        suptitle_y = 0.96
    elif len(facet_condition_list) == 3:
        suptitle_y = 0.98
    elif len(facet_condition_list) == 2:
        suptitle_y = 1.05
    elif len(facet_condition_list) == 1:
        suptitle_y = 1.1
    else: 
        suptitle_y = 0.96
    fig.suptitle(title if title != 'default' else f"{chip_id} {make_metric_name_look_nice(facet_by)} normalized coverage for all cells", y=suptitle_y, x = 0.5, size='large')

    legend_handles = []  # List to collect legend handles

    for i, condition in enumerate(facet_condition_list):
        df_subset = df_metrics[df_metrics[facet_by] == condition]

        #ax = axes[i] if len(facet_condition_list) > 1 else axes
        #lplot = sns.lineplot(data=df_subset, x='gc_content', y='normalized_coverage', color=pal[condition], units='cell_id', estimator=None, alpha=0.03, ax=ax)
        if len(facet_condition_list) == 1:
            lplot = sns.lineplot(data=df_subset, x='gc_content', y='normalized_coverage', color=pal[condition], units='cell_id', estimator=None, alpha=0.03)
            refplot = sns.lineplot(data=df_subset, x='gc_content', y='reference_gc_content', color = 'black')
        else:
            lplot = sns.lineplot(ax=axes[i], data=df_subset, x='gc_content', y='normalized_coverage', color=pal[condition], units='cell_id', estimator=None, alpha=0.03)
            refplot = sns.lineplot(ax=axes[i], data=df_subset, x='gc_content', y='reference_gc_content', color = 'black')
        lplot.set_ylim(0, 2)
        lplot.set_ylabel('Normalized coverage', fontsize=14)
        lplot.set_xlabel('GC content', fontsize=14)
        lplot.axhline(y=1.0, color='red')
        
        # Add handle for condition to the list
        condition_leg = mlines.Line2D([], [], color=pal[condition], markersize=15, label=condition)
        legend_handles.append(condition_leg)
        fig = lplot.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
    
    reference_leg = mlines.Line2D([], [], color='black', markersize=15, label='Reference GC Content')
    legend_handles.append(reference_leg)
    # Create a single legend with all the handles outside the plot on the upper right side
    plt.legend(handles=legend_handles, bbox_to_anchor=(1, 1))

    if save_png:
        plt.savefig(save_png, bbox_inches='tight', dpi=600)
        plt.close('all')
    if save_pdf:
        save_pdf.savefig(bbox_inches='tight')
        plt.close('all')


        
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
    return df_metrics




def make_wgs_plots(collated_input_dir, df, outdir, save_pdf = 'No', make_all = True):
    '''makes relevant plots to visually display wgs metrics results from a sequencing run - meant to be used after the bam_metrics pipeline
    collated_input_dir: directory where you have your collated output from the bam_metrics file
    df: dataframe you should get from running merge_sequencing_csvs_with_metrics_csvs()
    outdir: where you want your nice plot deposited'''
    chipid = df.chip_id.unique()[0]
    make_heatmap(df, 'mean_coverage_depth', outdir, plate_size=72, save_pdf=save_pdf)
    make_cellcondition_subplot_heatmaps(df, 'mean_coverage_depth', outdir, plate_size=72, save_pdf=save_pdf, change_scale = False, title = f"{chipid} coverage depth - fixed scale")
    make_cellcondition_subplot_heatmaps(df, 'mean_coverage_depth', outdir, plate_size=72, save_pdf=save_pdf, change_scale = True, title = f"{chipid} coverage depth - scale by condition")
    if make_all == True:
        df_wgs_hist= make_wgs_hist_df(collated_input_dir, df)
        make_density_lineplot(df_wgs_hist, outdir, "coverage", "high_quality_coverage_count", save_pdf=save_pdf, log = True)

              
            
def make_breadth_depth_subplots(chip_id, df, outdir, save_pdf = False):
    handles = make_handles(df, 'cell_condition', palette = colour_dict, include_n = True, include_median = False)
    legend_dict = {
        'handles': handles,
        'title': 'N cell condition',
        'bbox': (1.02, 1.0)}
    fig, axes = plt.subplots(1, 2, figsize = (10,5))
    plt.suptitle(f"{chip_id} coverage breadth and depth")
    swarm_subplot(df, 'cell_condition', 'mean_coverage_depth', axes, ax = 0,  palette = colour_dict, title = 'Depth')
    swarm_subplot(df, 'cell_condition', 'coverage_breadth', axes, ax = 1,  palette = colour_dict, title = 'Breadth', legend = legend_dict)
    plt.subplots_adjust(wspace=0.4, top=0.85)
    
    if save_pdf == False:
        plt.savefig(os.path.join(outdir, f'{chip_id}_depth_breadth_boxplot.png'), bbox_inches = 'tight')
        plt.close('all')
    else:
        save_pdf.savefig(bbox_inches = 'tight')
        plt.close('all')



def make_library_qc_plots(collated_input_dir, df, outdir):
    '''makes one pdf with key library metrics plots 
    collated_input_dir: directory with collated gatk metrics csvs produced by bam_metrics pipeline
    iso_csv: sequencing csv that came with isolatrix run 
    fastq_screen_csv: path to the collated fastq_screen csv produced by the fastq_metrics pipeline
    outdir: path to directory for plots to be created in'''
    if len(df.chip_id.unique()) == 1:
        chip_id = df.chip_id.unique()[0]
        outfile = os.path.join(outdir, f'{chip_id}_QC_plots.pdf')
        with PdfPages(outfile) as out_pdf:

            #make_breadth_depth_subplots(chip_id, df, outdir, save_pdf = out_pdf)
            make_wgs_plots(collated_input_dir, df, outdir, save_pdf = out_pdf, make_all = False)

            #make_mapping_plots(df, outdir, save_pdf = out_pdf, make_all = True)
            #make_insert_plots(collated_input_dir, df, outdir, save_pdf=out_pdf, make_all = False)
            #make_duplicate_plots(df, outdir, save_pdf = out_pdf)
            #make_gc_plots(collated_input_dir, df, outdir, save_pdf = out_pdf, make_all = False)
            #make_contam_boxplot(fastq_screen_csv, iso_csv, save_pdf= out_pdf, make_all = False)
        print('Successfully made key library metrics plots!') 
    else:
        print('QC plots can only be made for one chip at a time, and you have 2 chip ids in your dataframe')









def make_contam_df(fastq_screen_csv, iso_csv):
    fastq_screen_df = pd.read_csv(fastq_screen_csv)
    iso_df = pd.read_csv(iso_csv)
    iso_df['cell_id'] = (iso_df['sample_id'] + '-' + iso_df['chip_id'] + '-R' + iso_df['row'].apply(str) + '-C' + iso_df['column'].apply(str))
    merged = fastq_screen_df.merge(iso_df, how = 'left', on = 'cell_id')

    merged = merged.loc[merged['variable'].isin(['percent_one_hit_one_genome', 'percent_multiple_hits_one_genome', 'percent_hit_no_genomes'])]
    #merged['number_of_reads'] = merged['value']
    return merged



def make_contam_df(fastq_screen_csv, iso_csv):
    fastq_screen_df = pd.read_csv(fastq_screen_csv)
    iso_df = pd.read_csv(iso_csv)
    iso_df['cell_id'] = (iso_df['chip_id'] + '-' + iso_df['pool_id'] + '-R' + iso_df['row'].apply(str) + '-C' + iso_df['column'].apply(str))
    merged = fastq_screen_df.merge(iso_df, how = 'left', on = 'cell_id')

    merged = merged.loc[merged['variable'].isin(['percent_one_hit_one_genome', 'percent_multiple_hits_one_genome', 'percent_hit_no_genomes'])]
    #merged['number_of_reads'] = merged['value']
    return merged


def make_contam_boxplot(fastq_screen_csv, iso_csv, title = 'A130821A Livecells non-human contamination'):
    '''makes a boxplot displaying collated output from fastqscreen to indicate contamination in run
    fastq_screen_csv: str, is the path to the csv output by collate_fastq_screen function
    iso_csv: str, path to the csv given for the sequencing run '''
    df = make_contam_df(fastq_screen_csv, iso_csv)
    chip_id = df.chip_id.unique()[0]
    pal = {'percent_one_hit_one_genome': '#1f77b4', 'percent_multiple_hits_one_genome': '#ff7f0e', 'percent_hit_no_genomes': 'grey'}
    fig, ax = plt.subplots()
    plot = sns.barplot(data = df, x = 'genome', y = 'value', hue = 'variable', palette = pal, ax = ax)
    ax.set_title(title, size = 'medium')
    ax.set_ylabel('percent mapping')
    ax.set_xlabel('Genome')
    ax.get_legend().remove()
    plt.xticks(rotation = 90)
    
    one_hit_leg = mlines.Line2D([], [], color='#1f77b4',
                            markersize=15, label='One hit one genome')
    multi_hit_leg = mlines.Line2D([], [], color='#ff7f0e',
                    markersize=15, label='Multiple hits one genome')
    fig.legend(handles=[one_hit_leg, multi_hit_leg], bbox_to_anchor=(0.9,0.9))
    
    plt.savefig('fastq_screen_test.png', bbox_inches = 'tight')
plt.close('all')


def _contam_subplot(df, ax = 0, title = ' '):
    pal = {'percent_one_hit_one_genome': '#1f77b4', 'percent_multiple_hits_one_genome': '#ff7f0e', 'percent_hit_no_genomes': 'grey'}
    plot = sns.barplot(data = df, x = 'genome', y = 'value', hue = 'variable', palette = pal, ax = ax)
    ax.set_title(title, size = 'medium')
    ax.set_ylabel('percent mapping')
    ax.set_xlabel('Genome')
    ax.get_legend().remove()


def make_cellcondition_subplot_contamplots(fastq_screen_csv, iso_csv, outdir, title = 'default'):
    df = make_contam_df(fastq_screen_csv, iso_csv)
    chip_id = df.chip_id.unique()[0]
    if title == 'default': 
        title = f"{chip_id} read mapping"
    else:
        title = title 

    fig, axes = plt.subplots(2,2, figsize=(13,10))
    plt.suptitle(title, size= 'x-large', x = 0.3)
    i = 0
    for condition in df.cell_condition.unique():
        condition_subset = df[df["cell_condition"] == condition]
        subtitle = f"{condition}"
        _contam_subplot(condition_subset, ax = axes[i // 2, i % 2], title = subtitle)
        i += 1

    one_hit_leg = mlines.Line2D([], [], color='#1f77b4',
                            markersize=15, label='One hit one genome')
    multi_hit_leg = mlines.Line2D([], [], color='#ff7f0e',
                    markersize=15, label='Multiple hits one genome')
    fig.legend(handles=[one_hit_leg, multi_hit_leg], bbox_to_anchor=(0.9,1.05))
    plt.subplots_adjust(wspace=0.2, hspace=0.9, top=0.90)

    # Rotate x-axis tick labels for all subplots
    for ax in axes.flat:
        plt.sca(ax)
        plt.xticks(rotation=80, size = 'x-small')

    plt.savefig(f"{outdir}/{chip_id}_fastq_screen.png", bbox_inches = 'tight', dpi = 600)
    plt.close('all')



