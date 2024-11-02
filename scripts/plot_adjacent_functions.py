import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines
import pandas as pd

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



#print(sns.color_palette("Paired").as_hex())


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


def replace_condition_with_condition_n_in_df(df, df_condition, reverse = False):
    ''' Adds the number of cells that are included within a given condition
    df: dataframe
    df_condition: dataframe columns with condition you want to label with cells (ex. cell_condition)
    reverse: set to true if you want to remove n from a column where it has already been changed'''
    for condition in df[df_condition].unique():
        n = len(df[df[df_condition] == condition])
        if reverse == False:
            df.replace(condition, f'{condition} (n={n})', inplace = True)
        else:
            condition_no_n = condition.split(' (n=')[0]
            df.replace(condition, f'{condition_no_n}', inplace = True)



def update_pal_with_n(df, condition = 'cell_condition', palette = colour_dict):
    '''Use this to update the cell_condition palette if you have replaced the values of cell conditions (e.x. cell condition is now LiveCell (n = 21))'''
    new_pal = {}
    new_cell_conditions = df[condition].unique()
    for key in palette.keys():
        for new_condition in new_cell_conditions:
            if key in new_condition:
                new_pal[new_condition] = palette[key] 
    return new_pal




def make_handles(df, condition_metric, palette = False, include_n = True, include_median = False):
    '''Make new handles for a plot legend 
    df = dataframe 
    condition_metric: the condition you're splitting over and making the legend for (ex. chip_id)
    palette = the palette you're colouring your different conditions with - if false is seabron colours 
    include_n: Set to true if you want an (n=) included in the handles as a cell count 
    include median: this will need to be set for another column in your df (the dependant variable you're plotting - presumably on the y axis). Otherwise takes false'''
    handles = []
    i = 0
    condition = []
    for cond in df[condition_metric].unique():
        condition.append(cond)
    condition = sorted(condition)
    for cond in condition:
        if palette != False:
            colour = palette[cond]
        else:
            colour = sns.color_palette().as_hex()[i]
            
        if include_n == True and include_median != False:
            raise ValueError("You can't have the n of cells and the median of your variable in the legend. Please choose one or the other")
        elif include_n == True and include_median == False: 
            n = len(df[df[condition_metric] == cond])
            handle = mlines.Line2D([], [], color=colour,
                    markersize=25, label=str(f"{cond} (n={n})"))
        elif include_median != False and include_n == False:
            med = round((df[df[condition_metric] == cond])[include_median].median(), 5)
            handle = mlines.Line2D([], [], color=colour,
                    markersize=25, label=str(f"{cond} (Med={med})"))
        else:
            handle = mlines.Line2D([], [], color=colour,
                    markersize=25, label=str(f"{cond}"))
        handles.append(handle)
        i += 1
    return handles 



def make_handles_with_pal(pal):
    handles = []
    for condition in pal.keys():
        handle = plt.Circle((0.5, 0.5), 0.1, color=pal[condition], label=f"{make_metric_name_look_nice(condition)}")
        handles.append(handle)
    return handles 


def order_df_by_cell_condition(df):
    best_order = ['LiveCell', 'gDNA', 'NCC', 'NTC']
    existing_order = []
    for condition in best_order:
        if condition in df.cell_condition.unique():
            existing_order.append(condition)
    df['cell_condition'] = pd.Categorical(df['cell_condition'], categories=existing_order, ordered=True)
    df = df.sort_values(by='cell_condition')
    return df



