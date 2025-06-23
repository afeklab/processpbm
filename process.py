from utils import *

path1 = 'data/488nm_800_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
path2 = 'data/488nm_1000_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr'
path_fasta = 'data/086902_D_Fasta_20220417.txt'
f1 = 'F488 Median'
f2 = 'F488 Median'
path_result = None
ll = 200
lh = 40_000
radius = 7

# Set path_result directory tail as the one of path1 by default
if path_result == None:
    path_result = split(path1)[0]
path_masliner_png = join(path_result, 'masliner.png')
path_normalize_png = join(path_result, 'normalize.png')
path_processed_csv = join(path_result, 'processed.csv')
path_log = join(path_result, 'log.txt')
path_combinatorial = join(path_result, 'combinatorial.txt')

# Read GPRs
gpr1 = read_gpr(path1, f1) 
gpr2 = read_gpr(path2, f2) 

# Assert that the two GPRs are of the same layout and the same flags
assert np.all(gpr1[COLUMNS] == gpr2[COLUMNS])

# Perform masliner
y1, y2 = gpr1[f1].values, gpr2[f2].values
adj, slope, intercept, r, p, se, reg_num = masliner(y1, y2, ll, lh)
plot_masliner(y1, y2, adj, path_masliner_png)

# Create processed DataFrame with adjusted intensities
df = gpr2.drop(f2, axis=1)
df[f'Low - {f1}'] = y1
df[f'High - {f2}'] = y2
df['Adj'] = adj

# Perform normalization
df = normalize(df, path_fasta, 'Adj', radius=radius)
plot_normalize(df, 'Adj', path_normalize_png)
df.to_csv(path_processed_csv)

# Create combinatorial file
df_comb = df[['Sequence', 'Norm']]
df_comb = df_comb[~ df_comb['Norm'].isna()]
df_comb.sort_values('Norm', ascending=False, inplace=True)
df_comb.to_csv(path_combinatorial, header=False, index=False, sep='\t')

# Save result log
s = ''
s += '# Input information\n\n'
s += f'Low-scan GPR path = {path1}\n'
s += f'High-scan GPR path = {path2}\n'
s += f'f1 = {f1}\n'
s += f'f2 = {f2}\n'
s += f'll = {ll}\n'
s += f'lh = {lh}\n\n'
s += f'Fasta path = {path_fasta}\n'
s += f'Radius = {radius}\n\n'
s += '# Linear regression information\n\n'
s += f'Number of total probes = {len(y1)}\n'
s += f'Number of probes used for regression = {reg_num}\n'
s += f'Slope = {slope}\n'
s += f'Intercept = {intercept}\n'
s += f'Pearson correlation coefficient = {r}\n'
s += f'p-value = {p}\n'
s += f'Standrad error of the estimated slope = {se}\n\n'
s += '# Output information\n\n'
s += f'Processed path = {path_processed_csv}\n'
s += f'Combinatorial path = {path_combinatorial}\n'
s += f'Masliner figure path = {path_masliner_png}\n'
s += f'Normalize figure path = {path_normalize_png}\n'
s += f'Log path = {path_log}\n\n'

with open(path_log, 'w') as file:
    file.write(s)