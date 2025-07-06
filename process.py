from utils import *
import argparse

parser = argparse.ArgumentParser(prog='process', 
                                 description='Combining two PBM scans into one using liner regression (optional) and filtering out spatial correlation using a moving window',
                                )
parser.add_argument('-g1', '--gpr1', type=str, default=None, help='A path to a GPR file. If two GPR files are to be provided, this file is the low-intensity fluorescent scan one. Otherwise, the path directs to the single GPR file to be considered and the script will not perform Masliner.')
parser.add_argument('-g2', '--gpr2', type=str, default=None, help='A path to a GPR file. If two GPR files are to be provided, this file is the high-intensity fluorescent scan one. Otherwise, the path can directs to the single GPR file to be considered and the script will not perform Masliner.')
parser.add_argument('-f', '--fasta', type=str, default=None, help='A path to the fasta file. The normalization script considers only probes with a sequence provided in the fasta file and does not provide a normalized signal for any other (e.g., controls provided by the manufacturer). Moreover, the combinatorial file includes only probes with sequences provided in the Fasta file.')
parser.add_argument('-i', '--intensity_column', type=str, default=None, help='The fluorescent intensity column in the GPR files to be considered.')
parser.add_argument('-o', '--output', type=str, default=None, help='The path of output, in the form of <dir path>/<file prefix>.')
parser.add_argument('-ll' , '--low_bound', type=int, default=2_000, help='By default, 200. The linear lower bound for Masliner. This is the low threshold of intensity to consider for the linear regression.')
parser.add_argument('-lh', '--high_bound', type=int, default=50_000, help='By default, 50,000. The linear high bound for Masliner. This is the high threshold of intensity to consider for the linear regression.')
parser.add_argument('-r', '--radius', type=int, default=7, help='By default, 7. The radius of the moving window is used in normalization. Use 0 to skip normalization.')
args = parser.parse_args()

gpr1_path = args.gpr1
gpr2_path = args.gpr2
fasta_path = args.fasta
intensity_col = args.intensity_column
ll = args.low_bound
lh = args.high_bound
radius = args.radius
output = args.output


# Initiate log with input information
s = ''
s += '# Input information\n\n'
s += f'GPR1 path = {gpr1_path}\n'
s += f'GPR2 path = {gpr2_path}\n'
s += f'intensity_col = {intensity_col}\n'
s += f'll = {ll}\n'
s += f'lh = {lh}\n\n'
s += f'Fasta path = {fasta_path}\n'
s += f'Radius = {radius}\n\n'


if output == None:
    raise 'No output path' 
path_masliner_png = output + '_masliner.png'
path_normalize_png = output + '_normalize.png'
path_processed_csv = output + '_processed.csv'
path_log = output + '_log.txt'
path_combinatorial = output + '_combinatorial.txt'

is_masliner = gpr1_path != None and gpr2_path != None
is_norm = radius > 0

if gpr1_path == None and gpr2_path == None:
    raise 'Must have at least a single GPR file provided'

elif is_masliner:

    # Read GPRs
    gpr1 = read_gpr(gpr1_path, intensity_col) 
    gpr2 = read_gpr(gpr2_path, intensity_col) 

    # Assert that the two GPRs are of the same layout and the same flags
    if not np.all(gpr1[COLUMNS] == gpr2[COLUMNS]):
        raise 'Two GPRs are not of the same layout'

    # Perform masliner
    y1, y2 = gpr1[intensity_col].values, gpr2[intensity_col].values
    adj, slope, intercept, r, p, se, reg_num = masliner(y1, y2, ll, lh)
    plot_masliner(y1, y2, adj, path_masliner_png)

    # Create processed DataFrame with adjusted intensities
    df = gpr2.drop(intensity_col, axis=1)
    df[f'Low - {intensity_col}'] = y1
    df[f'High - {intensity_col}'] = y2
    df['Adj'] = adj
    col_to_norm = 'Adj'

    s += '# Linear regression information\n\n'
    s += f'Number of total probes = {len(y1)}\n'
    s += f'Number of probes used for regression = {reg_num}\n'
    s += f'Slope = {slope}\n'
    s += f'Intercept = {intercept}\n'
    s += f'Pearson correlation coefficient = {r}\n'
    s += f'p-value = {p}\n'
    s += f'Standrad error of the estimated slope = {se}\n\n'

else:
    if gpr1_path == None:
        gpr =read_gpr(gpr2_path, intensity_col)
    else:
        gpr =read_gpr(gpr1_path, intensity_col)
    df = gpr.copy()
    col_to_norm = intensity_col

# Merge df with fasta sequences (nan for rows with no correspondings in fasta)
df = pd.merge(df, read_fasta(fasta_path), how='left', on='ID')
 
if is_norm: 
    
    # Perform normalization
    final_col = 'Norm'
    df = normalize(df, col_to_norm, radius=radius)
    plot_normalize(df, col_to_norm, path_normalize_png)
    df.to_csv(path_processed_csv)
    mask_comb = ~ df[final_col].isna()

else:

    final_col = col_to_norm
    mask_flags, mask_seq = get_norm_masks(df)
    mask_comb = mask_flags * mask_seq

# Create combinatorial file
df_comb = df[['Sequence', final_col]][mask_comb]
df_comb.sort_values(final_col, ascending=False, inplace=True)
df_comb.to_csv(path_combinatorial, header=False, index=False, sep='\t')

s += '# Output information\n\n'
if is_masliner:
    s += f'Masliner figure path = {path_masliner_png}\n'
if is_norm:
    s += f'Normalize figure path = {path_normalize_png}\n'
s += f'Processed path = {path_processed_csv}\n'
s += f'Combinatorial path = {path_combinatorial}\n'
s += f'Log path = {path_log}\n\n'
with open(path_log, 'w') as file:
    file.write(s)
