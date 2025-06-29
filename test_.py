import subprocess
from os.path import join, expanduser # expand the tilde (~) character in a path string to the user's home directory. 
from utils import *
HOME = expanduser('~')
PBM_SUITE = 'PBM_analysis_suite_Sep2017'
PBM_SUITE_MASLINER = 'PBM_analysis_suite_Sep2017_masliner'
PBM_SUITE_NORMALIZE = 'PBM_analysis_suite_Sep2017_normalize'
f = 'F1 Median - B1'
ll = 200
lh = 40_000

def bash_unzip_suite(name):
    
    # Unzip PBM_analysis_suite
    # According to documentation: 
    # We specify to use the shell as the program to execute.
    # "If shell is True, it is recommended to pass args as a string rather than as a sequence."
    result = subprocess.run(f'unzip {PBM_SUITE}.zip -d {name}', shell=True, cwd=HOME, capture_output=True)
    subprocess.run(f'mv {name}/{PBM_SUITE}/* {name}/', shell=True, cwd=HOME)
    subprocess.run('mv masliner.txt masliner.pl', shell=True, cwd=join(HOME, name))

def test_maliner():

    bash_unzip_suite(PBM_SUITE_MASLINER)
    for version in [1, 2]:

        file1 = f'PBMscan_4x44k_v{version}_Alexa488_lp80pg80.gpr'
        file2 = f'PBMscan_4x44k_v{version}_Alexa488_lp100pg100.gpr'

        gpr1 = read_gpr(join(HOME, PBM_SUITE_MASLINER, 'test_files', file1), f)
        gpr2 = read_gpr(join(HOME, PBM_SUITE_MASLINER, 'test_files', file2), f)
        
        # Compute our masliner
        y1 = gpr1[f].values
        y2 = gpr2[f].values
        adj, _, _, _, _, _, _ = masliner(y1, y2, ll, lh)

        # Compute Perl masliner 
        perl_g1 = join('test_files', file1)
        perl_g2 = join('test_files', file2)
        subprocess.run('mkdir Temp_masliner', shell=True)
        perl_o = join(HOME, 'processpbm', 'Temp_masliner', 'masliner.gpr')

        result = subprocess.run(
            f'perl masliner.pl -g1 {perl_g1} -g2 {perl_g2} -o {perl_o} -ll {ll} -lh {lh} -f1 1 -f2 1',
            shell=True,
            cwd=join(HOME, PBM_SUITE_MASLINER),
            capture_output=True
        )
        perl_adj = read_gpr(join('Temp_masliner','masliner.gpr'), 'ADJBSI')['ADJBSI']

        # Assert the same result
        assert np.max(np.abs(perl_adj - adj)) == 0

        # Remove files
        subprocess.run('rm -rf Temp_masliner', shell=True)
    subprocess.run(f'rm -rf {join(HOME, PBM_SUITE_MASLINER)}', shell=True)

def test_normalize():


    bash_unzip_suite(PBM_SUITE_NORMALIZE)

    for version in [1, 2]:
    
        subprocess.run('mkdir Temp_normalize', shell=True)

        # Create Fasta file from probe_sequences in PBM suite
        fasta_file = f'4x44k_v{version}_sequences.txt'
        df_fasta = pd.read_csv(join(HOME, PBM_SUITE_NORMALIZE, 'probe_sequences', fasta_file), sep='\t', index_col=False)
        ID = np.array([f'>{ID_i}\n' for ID_i in df_fasta['ID']])
        seq = np.array([f'{seq_i}\n' for seq_i in df_fasta['Sequence']])
        ID_seq = np.char.add(ID, seq)
        ID_seq = ID_seq[(np.char.find(ID_seq, 'dBr') > -1) + (np.char.find(ID_seq, 'Ctrl') > -1)] # Mask only dBr sequences
        fasta_str = ''.join(ID_seq)
        with open(join('Temp_normalize', f'fasta_v{version}.txt'), 'w') as file:
            file.write(fasta_str)


        file = f'PBMscan_4x44k_v{version}_Alexa488_MaslinerOutput.gpr'
        gpr = read_gpr(join(HOME, PBM_SUITE_NORMALIZE, 'test_files', file), 'ADJBSI')

        
        gpr = pd.merge(gpr, read_fasta(join('Temp_normalize', f'fasta_v{version}.txt')), how='left', on='ID')
        df_norm = normalize(gpr, 'ADJBSI', corners=[20, 245, 150, 416])
        # plot_normalize(df_norm, 'ADJBSI')

        perl_i = join('test_files', file)
        perl_s = join('probe_sequences', fasta_file)
        perl_o = join(HOME, 'processpbm', 'Temp_normalize', 'normalize')
        result = subprocess.run(
            f'perl normalize_agilent_array.pl -i {perl_i} -s {perl_s} -o {perl_o}',
            shell=True,
            cwd=join(HOME, PBM_SUITE_NORMALIZE),
            capture_output=True
        )
        df_perl_norm = pd.read_csv(join('Temp_normalize', 'normalize_alldata.txt'), sep='\t', index_col=None)
        df_perl_norm['ID'] = df_perl_norm['ID'].astype(str)

        for c in ['Column', 'Row', 'ID']:
            assert np.all(df_norm[c] == df_perl_norm[c])

        assert np.max(np.abs(df_norm['Norm'] - df_perl_norm['Alexa488Adjusted'])) < 1e-9
        # plt.scatter(df_norm['Norm'], df_perl_norm['Alexa488Adjusted'])
        # plt.show()

        # Remove files
        subprocess.run('rm -rf Temp_normalize', shell=True)
    subprocess.run(f'rm -rf {join(HOME, PBM_SUITE_NORMALIZE)}', shell=True)