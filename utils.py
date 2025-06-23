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

def masliner(path1, path2, f1, f2, ll=200, lh=40_000, path_result=None):

    # Set path_result directory tail as the one of path1 by default
    if path_result == None:
        path_result = split(path1)[0]
    path_csv = join(path_result, 'masliner.csv')
    path_png = join(path_result, 'masliner.png')
    path_txt = join(path_result, 'masliner.txt')

    # Read GPRs
    gpr1 = read_gpr(path1, f1) 
    gpr2 = read_gpr(path2, f2) 

    # Assert that the two GPRs are of the same layout and the same flags
    assert np.all(gpr1[COLUMNS] == gpr2[COLUMNS])

    # Get fluorescent intensity
    y1, y2 = gpr1[f1].values, gpr2[f2].values

    # Mask range for linear regression and perform linear regression (PerformStraightRegression in original)
    mask = (y1 >= ll) * (y2 >= ll) * (y1 <= lh) * (y2 <= lh)
    slope, intercept, r, p, se = linregress(y1[mask], y2[mask])
    linreg_pred = (y1 * slope + intercept + 0.5).astype(int) # Round to integers (PerformAdjustment in original)
    adj = np.where(y2 > lh, linreg_pred, y2)
    gpr2.drop(f2, axis=1, inplace=True)
    gpr2[f'Low - {f1}'] = y1
    gpr2[f'High - {f2}'] = y2
    gpr2['Adj'] = adj

    # Save adjusted gpr file
    gpr2.to_csv(path_csv)

    # Plot result
    fig, axis = plt.subplots()
    axis.set_title('Masliner adjusted result')
    axis.scatter(y1, y2, label='original')
    axis.scatter(y1, adj, label='adjusted')
    axis.set_xlabel('Low-scan fluorescent intensity')
    axis.set_ylabel('High-scan fluorescent intensity')
    axis.legend()
    plt.savefig(path_png)

    # Save result log
    s = ''
    s += '# Masliner result\n\n'
    s += '## Input information\n\n'
    s += f'Low-scan GPR path = {path1}\n'
    s += f'High-scan GPR path = {path2}\n'
    s += f'f1 = {f1}\n'
    s += f'f2 = {f2}\n'
    s += f'll = {ll}\n'
    s += f'lh = {lh}\n\n'
    s += '## Linear regression information\n\n'
    s += f'Number of total probes = {len(y1)}\n'
    s += f'Number of probes used for regression = {np.sum(mask)}\n'
    s += f'Slope = {slope}\n'
    s += f'Intercept = {intercept}\n'
    s += f'Pearson correlation coefficient = {r}\n'
    s += f'p-value = {p}\n'
    s += f'Standrad error of the estimated slope = {se}\n\n'
    s += '## Output information\n\n'
    s += f'Masliner result path = {path_csv}\n'
    s += f'Masliner figure path = {path_png}\n'
    s += f'Masliner log path = {path_txt}\n\n'
    with open(path_txt, 'w') as file:
        file.write(s)

# Custom formatter: convert to thousands with 'k'
def k_formatter(x, pos):
    return f'{int(x/1000)}k'

def normalize(path, path_fasta, f='Adj', radius=7, custom_mask=True, path_result=None):

    # Set path_result directory tail as the one of path1 by default
    if path_result == None:
        path_result = split(path1)[0]
    path_csv = join(path_result, 'normalize.csv')
    path_comb = join(path_result, 'combinatorial.txt')
    path_png = join(path_result, 'normalize.png')
    path_txt = join(path_result, 'normalize.txt')

    # Get data
    if path[-4:] == '.gpr':
        df = read_gpr(path, f)
    else:
        df = pd.read_csv(path, index_col=0)
    
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

    df['Top'] = top_adj
    df['Bottom'] = bottom_adj
    df['Left'] = left_adj
    df['Right'] = right_adj

    # Normalize
    median_global = np.median(y)
    median_local = np.empty(len(df))
    size = np.empty(len(df))
    '''
    print('Running window')
    for i in tqdm(range(len(df))):
    
        t, b, l, r = df[['Top', 'Bottom', 'Left', 'Right']].iloc[i].values
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
    df['Local window median'] = median_local
    df['Window size'] = size
    df['Norm'] = norm
    df.to_csv(path_csv)
    '''
    df = pd.read_csv(path_csv)

    df_comb = df[['Sequence', 'Norm']]
    df_comb = df_comb[~ df_comb['Norm'].isna()]
    df_comb.sort_values('Norm', ascending=False, inplace=True)
    df_comb.to_csv(path_comb, header=False, index=False, sep='\t')

    # Plots result
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    arrays = [y_arr, df['Norm'].values.reshape(row_max, col_max), mask_arr]
    titles = [f, 'Norm', 'Norm mask (binary)']

    for i, ax, arr, title in zip(range(3), axes, arrays, titles):
        im = ax.matshow(arr)
        ax.set_title(title)

        if i < 2:
            
            # Create a divider for each axis
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(im, cax=cax)
            cbar.ax.yaxis.set_major_formatter(FuncFormatter(k_formatter))

    plt.tight_layout()
    plt.savefig(path_png)

    # Save result log
    s = ''
    s += '# Normalize result\n\n'
    s += '## Input information\n\n'
    s += f'Data path = {path}\n'
    s += f'Fasta path = {path_fasta}\n'
    s += f'f = {f}\n'
    s += f'Radius = {radius}\n\n'
    s += '## Output information\n\n'
    s += f'Normalize result path = {path_csv}\n'
    s += f'Combinatorial path = {path_comb}\n'
    s += f'Normalize figure path = {path_png}\n'
    s += f'Normalize log path = {path_txt}\n\n'
    with open(path_txt, 'w') as file:
        file.write(s)
    
path1 = 'data/488nm_800_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
path2 = 'data/488nm_1000_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
# masliner(path1, path2, 'F488 Median', 'F488 Median')
normalize('data/masliner.csv', 'data/086902_D_Fasta_20220417.txt')