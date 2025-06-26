import subprocess
from os.path import join, expanduser # expand the tilde (~) character in a path string to the user's home directory. 
from utils import *
HOME = expanduser('~')
PBM_SUITE = 'PBM_analysis_suite_Sep2017'
PBM_SUITE_MASLINER = 'PBM_analysis_suite_Sep2017_masliner'

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
    f = 'F1 Median - B1'
    ll = 200
    lh = 40_000

    for ver in [1, 2]:

        # Read GPRs
        gpr1_file = f'PBMscan_4x44k_v{ver}_Alexa488_lp80pg80.gpr'
        gpr2_file = f'PBMscan_4x44k_v{ver}_Alexa488_lp100pg100.gpr'

        gpr1 = read_gpr(join(HOME, PBM_SUITE_MASLINER, 'test_files', gpr1_file), f)
        gpr2 = read_gpr(join(HOME, PBM_SUITE_MASLINER, 'test_files', gpr2_file), f)


        # Compute our masliner
        y1 = gpr1[f].values
        y2 = gpr2[f].values
        adj, _, _, _, _, _, _ = masliner(y1, y2, ll, lh)

        # Compute Perl masliner 
        perl_g1 = join('test_files', gpr1_file)
        perl_g2 = join('test_files', gpr2_file)
        subprocess.run('mkdir Temp_masliner', shell=True)
        perl_o = join(HOME, 'processpbm', 'Temp_masliner', 'masliner.gpr')

        result = subprocess.run(
            f'perl masliner.pl -g1 {perl_g1} -g2 {perl_g2} -o {perl_o} -ll {ll} -lh {lh} -f1 1 -f2 1',
            shell=True,
            cwd=join(HOME, PBM_SUITE_MASLINER),
            capture_output=True
        )
        perl_adj = read_gpr('Temp_masliner/masliner.gpr', 'ADJBSI')['ADJBSI']

        # Assert the same result
        assert np.max(np.abs(perl_adj - adj)) == 0

        # Remove files
        subprocess.run('rm -rf Temp_masliner', shell=True)
    subprocess.run(f'rm -rf {join(HOME, PBM_SUITE_MASLINER)}', shell=True)