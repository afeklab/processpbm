import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from os.path import split, join

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

def normalize(path, path_fasta, f='Adj', radius=7):

    if path[-4:] == '.gpr':
        df = read_gpr(path, f)
    else:
        df = pd.read_csv(path, index_col=0)
    
    # Merge df with fasta sequences (nan for rows with no correspondings in fasta)
    fasta = read_fasta(path_fasta)
    df = pd.merge(df, fasta, how='left', on='ID')

    # Set boundaries for moving window
    df['Top'] = df['Row'] - radius
    df['Bottom'] = df['Row'] + radius
    df['Left'] = df['Column'] - radius
    df['Right'] = df['Column'] + radius

    # Adjust 'out-of-boundaris' windows
    # mask_no_seq = df['Sequence'].isna()
    print(df.iloc[:20])

path1 = 'data/488nm_800_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
path2 = 'data/488nm_1000_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
# masliner(path1, path2, 'F488 Median', 'F488 Median')
normalize('data/masliner.csv', 'data/086902_D_Fasta_20220417.txt')