import matplotlib.pyplot as plt
import seaborn as sns
import sys
from plot_adjacent_functions import *
sys.path.insert(0, '/projects/steiflab/research/hmacdonald/other_code/comparisons')
from compare_library_functions import *

#set colour palette and theme for plots
pal = sns.color_palette()
sns.set_palette(pal)
sns.set_context("talk")


def make_swarm_boxplot(df, xmetric, ymetric, hue = None, title = 'default', log = False, legend = True, include_n_in_legend = True, handles = 'default', bbox = (1,1), outplot = 'swarm_boxplot.png', palette = 'default', xlabel = ' ', ylabel = 'default', x_rotate = False, include_n_in_x = False, extra_code=None, savefig = True, dpi = 100):
    if hue == None:
        hue = xmetric
    if palette == 'default':
        palette = make_palette(df, xmetric)
    categories = df[xmetric].unique()
    if include_n_in_x == True:
        for category in categories:
            len_df_category = len(df[df[xmetric] == category])
            new_df = df.replace(category, f'{category} (n = {len_df_category})')
    else:
        new_df = df
    if include_n_in_x == True and include_n_in_legend == True:
        palette = update_pal_with_n(new_df, condition = xmetric, palette = palette)

    if include_n_in_x == False and include_n_in_legend == True:
        handles = make_handles(new_df, xmetric, palette = palette, include_n = True, include_median = False)

    fig, ax = plt.subplots()
    if log == True:
        new_df[ymetric] = new_df[ymetric] + 0.01
        ax.set_yscale('log')
    vi = sns.boxplot(data=new_df, y=ymetric, x = xmetric,  color = 'lightgrey')
    if legend == True:
        sw = sns.stripplot(data=new_df, y=ymetric, x = xmetric, hue = hue, palette = palette, size = 3, edgecolor = 'black')
        
        if handles == 'default':
            sw.legend(bbox_to_anchor=bbox)
        else:
            sw.legend(handles = handles, bbox_to_anchor=bbox)
    else:
        sw = sns.stripplot(data=new_df, y=ymetric, x = xmetric, hue = hue, palette = palette, size = 3, edgecolor = 'black', legend = False)
     
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
    if savefig == True:
        plt.savefig(outplot,  bbox_inches = 'tight', dpi = dpi)
        plt.close('all')


def reset_xmetric_names(df, xmetric):
    """Resets a metric if n= or med = has been added to it - helpful if you're making multiple swarm boxplots in a row"""
    for cat in df[xmetric].unique():
        new_cat = cat.split(' ')[0]
        df.replace(cat, new_cat, inplace = True)

        
def swarm_subplot(df, xmetric, ymetric, axes, ax = 0, title = ' ', ylabel = ' ', xlabel = ' ', log = False, palette = False, legend = 'None', x_rotate = False):
    '''Make a swarm boxplot in one panel of a faceted plot
    ax: the axis you want to plot into 
    axes: the axes object created by calling fig,axes = plt.subplots() - this needs to be called before running this function
    If you want a legend, feed legend a dictionary in this format:
    {handles: handles,
    title: 'title for legend'
    bbox: (1,1)}'''

    sns.boxplot(data=df, x=xmetric, y=ymetric, color = 'lightgray', ax=axes[ax])
    if palette != False:
        sns.stripplot(data=df, x=xmetric, y=ymetric, hue = xmetric, dodge = False, size = 3, edgecolor = 'black', ax=axes[ax], palette = palette)
    else:
        sns.stripplot(data=df, x=xmetric, y=ymetric, hue = xmetric, dodge = False, size = 3, edgecolor = 'black', ax=axes[ax])

    if legend != 'None':
        axes[ax].legend(handles = legend['handles'], title = legend['title'], bbox_to_anchor= legend['bbox'])
    else:
        axes[ax].get_legend().set_visible(False)

    axes[ax].set_title(title)
    if log == True:
        axes[ax].set_yscale('log')
    axes[ax].set_ylabel(ylabel)
    axes[ax].set_xlabel(xlabel)
    if x_rotate == True:
        axes[ax].tick_params(axis='x', rotation=45)  # Rotate x-axis labels


    

def make_density_lineplot(dataframe, outdir, xmetric, ymetric, colour_by = 'cell_condition', palette = colour_dict, split_by = 'cell_id', save_pdf=False, outplot = 'density_lineplot.png', log = False, bbox = (1,1)):
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
            
        if save_pdf == False:
            plt.savefig(outplot, dpi = 600, bbox_inches = 'tight')
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
            
            
#def make_subplots():


