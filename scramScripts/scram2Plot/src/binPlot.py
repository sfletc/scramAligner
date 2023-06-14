import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scram2Plot import *

def get_bin_data(srp, bin_size):
    bin_no = math.ceil(srp.ref_len/bin_size)
    
    # bins is a lists of lists of len bin no., with a f and r list in each list
    bins = np.array([[[0.]*srp.replicates, [0.]*srp.replicates] for _ in range(bin_no)])
    for i in srp.all_alignments:
        bin_count = int(i.position/bin_size)
        if i.strand == "+":
            bins[bin_count][0]+=i.indv_alignments
        else:
            bins[bin_count][1]+=i.indv_alignments
    return bins

def calculate_means_and_ses(bins):
    means = np.mean(bins, axis=2)
    ses = np.std(bins, axis=2) / np.sqrt(bins.shape[2])
    return means, ses

def create_dataframe(bins, means, ses):
    data = []
    for i in range(bins.shape[0]):
        data.append([i, 'first', means[i, 0], ses[i, 0]])
        data.append([i, 'second', -means[i, 1], ses[i, 1]])
    df = pd.DataFrame(data, columns=['bin', 'list', 'mean', 'se'])
    return df

def plot_standard_errors(df, color='gray', alpha=0.3):
    for list_name in ['first', 'second']:
        plt.fill_between(df[df['list'] == list_name]['bin'], 
                         (df[df['list'] == list_name]['mean'] - df[df['list'] == list_name]['se']).values, 
                         (df[df['list'] == list_name]['mean'] + df[df['list'] == list_name]['se']).values, 
                         color=color, alpha=alpha)
    return df

def create_plot(bins, color='gray', alpha=0.7):
    means, ses = calculate_means_and_ses(bins)
    df = create_dataframe(bins, means, ses)
    df = plot_standard_errors(df, color, alpha)
    return df

def show_plot(dataframes, scale_factor=1, filename=None, y_min=None, y_max=None, cols={}):
    plt.axhline(0, color='black')
    if y_min is None or y_max is None:  # If y-axis limits are not provided
        y_abs_max = max(max(df['mean'] + df['se']) for df in dataframes) * 1.1  # Increased limits by 10%
        if y_min is None:
            y_min = -y_abs_max
        if y_max is None:
            y_max = y_abs_max
    plt.ylim(y_min, y_max)
    num_bins = max(df['bin'].max() for df in dataframes) + 1
    ticks = np.linspace(0, num_bins - 1, 10).astype(int)
    scaled_ticks = ticks * scale_factor
    formatted_ticks = [format(int(tick), ',') for tick in scaled_ticks]  # Format with commas
    plt.xticks(ticks, formatted_ticks, rotation=45)  # Rotate labels 45 degrees
    plt.xlabel("Reference position (bp)")  # Label for x-axis
    plt.ylabel("Reads per million reads")  # Label for y-axis
    create_custom_legend(cols)  # create custom legend
    if filename:
        plt.savefig(filename, bbox_inches='tight')  # Adjust bounding box to fit rotated labels
    else:
        plt.show()

def create_custom_legend(colors_dict, alpha=0.7):
    # Order by converting keys to integers and sorting
    colors_ordered = {k: colors_dict[k] for k in sorted(colors_dict, key=int)}
    patches = [mpatches.Patch(color=color, label=str(label), alpha=alpha) for label, color in colors_ordered.items()]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left')
    
def plot_bin_plot(srps, lens, bin_size, filename=None, y_min=None, y_max=None):
    plt.figure(figsize=(10,6))
    cols={24:'darkgreen', 21:'#984ea3', 22:'#3374FF',
                  18:'#e41a1c', 19:'#33D4FF', 20:'#4daf4a',
                  23:'#ffff33', 25:'#a65628', 26:'#f781bf'}
    dfs = []
    pos = 0
    for i in srps:
        dfs.append(create_plot(get_bin_data(i,bin_size),cols[lens[pos]]))
        pos+=1
    show_plot(dfs, bin_size, filename, y_min, y_max, cols)

# Load single_ref_profiles from scram2Plot
def load_ref_profiles(file_prefix, lens):
    rp=[]
    for i in lens:
        file_name = "{}_{}.csv".format(file_prefix, i)
        rprofile=RefProfiles()
        rprofile.load_single_ref_profiles(file_name)
        rp.append(rprofile)
    return rp

def bin_plot(file_prefix, lens, headers, bin_size, out_prefix=None, y_min=None, y_max=None):
    rps = load_ref_profiles(file_prefix, lens)
    
    for header in headers:
        print(header)
        srps = []
        for i in rps:
            try:
                srp = i.single_ref_profiles[header]
                srps.append(srp)
            except:
                pass
        if out_prefix is not None:
            filename = "{}_{}.png".format(out_prefix, header.split()[0])
        plot_bin_plot(srps, lens, bin_size, filename, y_min, y_max)



