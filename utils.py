import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from tqdm import tqdm
from os.path import split, join
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter

COLUMNS = ['Block', 'Column', 'Row', 'Name', 'ID', 'Flags']

def read_gpr(path, f):

    
    # Get the index of the header of the GPR without reading the entire file    
    with open(path, 'r') as file:

        # Read line by line until reaching COLUMNS and f (only one line is in memory at a iteration)     
        header_i = 0
        for line in file:
            
            # If all columns in the line
            if np.all([c in line for c in COLUMNS + [f]]):
                break
            else: 
                header_i += 1
    
    # Read and return GPR DataFrames using pandas and the header index
    gpr = pd.read_csv(path, sep='\t', header=header_i)[COLUMNS + [f]]
    gpr['Name'] = gpr['Name'].astype(str)
    gpr['ID'] = gpr['ID'].astype(str)
    return gpr

def read_fasta(path):
    fasta = pd.read_csv(path, lineterminator='>', header=None)[0] # Select the single colum to get a Series
    fasta = fasta.str.extract(r'(?P<ID>.*)\n(?P<Sequence>.*)\n') # Extract ID and sequence using regular expressions (see function documentation for explanation)
    return fasta

def masliner(y1, y2, ll=200, lh=40_000):

    # Mask range for linear regression and perform linear regression (PerformStraightRegression in original)
    mask = (y1 >= ll) * (y2 >= ll) * (y1 <= lh) * (y2 <= lh)
    reg_num = np.sum(mask)
    slope, intercept, r, p, se = linregress(y1[mask], y2[mask])
    linreg_pred = (y1 * slope + intercept + 0.5).astype(int) # Round to integers (PerformAdjustment in original)
    adj = np.where(y2 > lh, linreg_pred, y2)
    return adj, slope, intercept, r, p, se, reg_num

def normalize(df, path_fasta, f='Adj', radius=7, custom_mask=True):

    column = df['Column'].values
    row = df['Row'].values
    y = df[f].values
    mask_flags = (df['Flags'] > -100).values    
    col_max, row_max = np.max(column), np.max(row)

    # Check reshaped column is in the tile format of [[1, 2, 3, ...], [1, 2, 3, ...], ...] = np.tile(np.arange(col_max) + 1, (row_max, 1))
    assert np.all(column.reshape(row_max, col_max) == np.tile(np.arange(col_max) + 1, (row_max, 1)))

    # Check reshaped row is in the repeated format of [[1, 1, 1, ...], [2, 2, 2, ...], ...] = np.tile(np.arange(row_max) + 1, (col_max, 1)).T (similar to the column format but transposed)
    assert np.all(row.reshape(row_max, col_max) == np.tile(np.arange(row_max) + 1, (col_max, 1)).T)   
    
    # Merge df with fasta sequences (nan for rows with no correspondings in fasta)
    fasta = read_fasta(path_fasta)
    df = pd.merge(df, fasta, how='left', on='ID')
    mask_seq = ~ df['Sequence'].isna().values
    
    # Set intensity and mask in a 2-dimensional format of the array
    y_arr = y.reshape(row_max, col_max)
    mask_arr = (mask_seq * mask_flags * custom_mask).reshape(row_max, col_max)
    
    # Set boundaries for moving window
    top = row - radius
    bottom = row + radius
    left = column - radius
    right = column + radius

    top_adj = top.copy()
    bottom_adj = bottom.copy()
    left_adj = left.copy()
    right_adj = right.copy()

    # Adjust deviation from array top
    con = top < 1
    top_adj = np.where(con, 1, top_adj)
    bottom_adj = np.where(con, 1 + 2 * radius, bottom_adj)

    # Adjust deviation from array bottom
    con = bottom > row_max
    top_adj = np.where(con, row_max - 2 * radius, top_adj)
    bottom_adj = np.where(con, row_max, bottom_adj)

    # Adjust deviation from array left
    con = left < 1
    left_adj = np.where(con, 1, left_adj)
    right_adj = np.where(con, 1 + 2 * radius, right_adj)

    # Adjust deviation from array right
    con = right > col_max
    left_adj = np.where(con, col_max - 2 * radius, left_adj)
    right_adj = np.where(con, col_max, right_adj)

    # Normalize
    median_global = np.median(y)
    median_local = np.empty(len(df))
    size = np.empty(len(df))

    for i, t, b, l, r in zip(range(len(df)), top_adj, bottom_adj, left_adj, right_adj):# Note that running an index, take row from the dataframe, convert it to nu,py values at each iteration is slow
        
        y_window = y_arr[t - 1: b, l - 1: r].flatten()
        mask_window = mask_arr[t - 1: b, l - 1: r].flatten()
        median_local[i] = np.median(y_window[mask_window])
        size[i] = np.sum(mask_window)

    # If a probe has a window greater than the half of the maximum, has a sequence and appropriate flag, use normalized value, otherwise, use original
    mask_norm = (size > (1 + 2 * radius) ** 2 / 2) * (mask_seq) * (mask_flags)
    norm = np.where(mask_norm, y * median_global / median_local, y)
    
    # If probe is flagged or with no sequence, assign nan
    norm = np.where((~ mask_flags) + (~ mask_seq), np.nan, norm)

    # Save result
    df['Top'] = top_adj
    df['Bottom'] = bottom_adj
    df['Left'] = left_adj
    df['Right'] = right_adj
    df['Local window median'] = median_local
    df['Window size'] = size
    df['Norm'] = norm

    return df
def plot_masliner(y1, y2, adj, path):
    
    fig, axis = plt.subplots()
    axis.set_title('Masliner adjusted result')
    axis.scatter(y1, y2, label='original')
    axis.scatter(y1, adj, label='adjusted')
    axis.set_xlabel('Low-scan fluorescent intensity')
    axis.set_ylabel('High-scan fluorescent intensity')
    axis.legend()
    plt.savefig(path)


def plot_normalize(df, f, path):
    
    column = df['Column'].values
    row = df['Row'].values
    col_max, row_max = np.max(column), np.max(row)

    y_arr = df[f].values.reshape(row_max, col_max)
    norm_arr = df['Norm'].values.reshape(row_max, col_max)

    k_formatter = lambda x, pos : f'{int(x/1000)}k' # Custom formatter: convert to thousands with 'k'
    
    # Plots result
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    arrays = [y_arr, norm_arr]
    titles = [f, 'Normalized']

    for ax, arr, title in zip(axes, arrays, titles):
        im = ax.matshow(arr)
        ax.set_title(title)

        # Create a divider for each axis
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(k_formatter))

    plt.tight_layout()
    plt.savefig(path)



