import numpy as np
import pandas as pd
from scipy.stats import linregress
from os.path import split, join

COLUMNS = ['Block', 'Column', 'Row', 'Name', 'ID', 'Flags']

def read_gpr(path, f):

    
    # Get the index of the header of the GPR without reading the entire file    
    with open(path, 'r') as file:

        # Read line by line until reaching the column y (only one line is in memory at a iteration)     
        header_i = 0
        for line in file:
            
            # If all columns in the line
            if np.all([c in line for c in COLUMNS + [f]]):
                break
            else: 
                header_i += 1
    
    # Read and return GPR using pandas and the header index
    gpr = pd.read_csv(path, sep='\t', header=header_i)[COLUMNS + [f]]
    gpr['Name'] = gpr['Name'].astype(str)
    gpr['ID'] = gpr['ID'].astype(str)
    return gpr

def masliner(path1, path2, f1, f2, ll=200, lh=40_000, path_result=None):

    # Set path_result directory tail as the one of path1 by default
    if path_result == None:
        path_result = split(path1)[0]

    # Read GPRs
    gpr1 = read_gpr(path1, f1) 
    gpr2 = read_gpr(path2, f2) 

    # Assert that the two GPRs are of the same layout and the same flags
    assert np.all(gpr1[COLUMNS] == gpr2[COLUMNS])

    # Get fluorescent intensity
    y1, y2 = gpr1[f1].values, gpr2[f2].values

    # Mask range for linear regression and perform linear regression (PerformStraightRegression)
    mask = (y1 >= ll) * (y2 >= ll) * (y1 <= lh) * (y2 <= lh)
    slope, intercept, r, p, se = linregress(y1[mask], y2[mask])
    linreg_pred = (y1 * slope + intercept + 0.5).astype(int) # Rounding to integers, see PerformAdjustment() in original

    adj = np.where(y2 > lh, linreg_pred, y2)
    gpr2.drop(f2, axis=1, inplace=True)
    gpr2[f'low: {f1}'] = y1
    gpr2[f'high: {f2}'] = y2
    gpr2['adj'] = adj

    # Save adjusted gpr file
    gpr2.to_csv(join(path_result, 'masliner.csv'))

path1 = 'data/488nm_800_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
path2 = 'data/488nm_1000_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'

masliner(path1, path2, 'F488 Median', 'F488 Median')