# Processing Protein Binding Microarray Data

## Introduction

The fluorescent signal of a Protein Binding Microarray (PBM) carries information about the binding specificities of a transcription factor (TF). 
However, the signal does not explicitly expose the binding specificity.
Often, experimental sources in the assay alter or truncate the signal.
The Python scripts here aim to remove such experimental variations, leaving the specific-binding signal from a PBM.
We designed the scripts to perform the same functionalities as provided originally in `PBM_analysis_suite_Sep2017` Perl scripts (Berger et al., 2009), which can be downloaded from the Bulyk lab website.
We reimplemented the Masliner and spatial-normalization functions, removed the features we do not usually perform, and made the code easier to read/utilize.


## Quick Start

1. Make sure your Python environment contains the installed packages of Numpy, Pandas, Matplotlib and Scipy.

2. Change directory in the terminal to the one that contains the `process.py`, e.g., in bash:
```
cd processpbm
```

3. Run:
```
python process.py [-g1 GPR1] [-g2 GPR2] [-f FASTA] [-i INTENSITY_COLUMN] [-o OUTPUT] [-ll LOW_BOUND] [-lh HIGH_BOUND] [-r RADIUS] 

```
### Arguments

`-g1` : str, optional if `g2` is provided

A path to a GPR file.
If two GPR files are to be provided, this file is the low-intensity fluorescent scan one.
Otherwise, the path directs to the single GPR file to be considered and the script will nor perform Masliner.

`-g2` : str, optional if `g1` is provided

A path to a GPR file.
If two GPR files are designed to be provided, this file is the high-intensity fluorescent scan one.
Otherwise, the path can directs to the single GPR file to be considered and the script will nor perform Masliner.

`-f` :  str

A path to the fasta file.
Fasta file is a text file that contains all probe IDs and sequences of the PBM design in fasta form (see example in `data`).
The normalization script considers only probes with sequence provided in the fasta file and does not provide normalized signal for any other (e.g., controls provided by manufacture).
Moreover, the combinatorial file includes only probes with sequences provided in the fasta file.
Agilent provides the fasta file in the appropriate format for all designs.

`-i` : str

The fluorescent intensity column in the GPR files should be considered.
Note that in the original scripts it is the median intensity substructed by the background.

`-o` : str

The path of output is in the form of `<dir path>/<file prefix>`.


`-ll` : int, optional

By default, 200.
The linear low bound for Masliner.
This is the low threshold of intensity to consider for the linear regression.

`-lh` : int, optional

By default, 40,000.
The linear high bound for Masliner.
This is the high threshold of intensity to consider for the linear regression.

`-r` : int, optional

By default, 7.
The radius of the moving window is used in normalization.
Use 0 to skip normalization.

### Output

TThe script saves multiple files with the directory and prefix file name provided in the output path argument.
The following are the file suffixes and their descriptions:

`_log.txt` : 
The log text file of the script run.

`_processed.csv` : 
Contains the following columns:
* Block (from original GPR)
* Column (from original GPR)
* Row (from original GPR)
* Name (from original GPR)
* ID (from original GPR)
* Flags (from original GPR)
* \<i\> (if Masliner was not performed) - fluorescent intensity column from original GPR.
* Low $^*$ - low fluorescent intensity column from original GPR.
* High $^*$ - high fluorescent intensity column from original GPR.
* Adj $^*$ - Masliner adjusted values.
* Sequence (from fasta file)
* Top $^†$ - The top row boundary of the moving window.
* Bottom $^†$ - The bottom row boundary of the moving window.
* Left $^†$ - The left column boundary of the moving window.
* Right $^†$ - The right column boundary of the moving window.
* Local window median $^†$ - The considered local median intensity of the window.
* Window size $^†$ - The number of probes in the window used for local median intensity calculation.
* Norm $^†$ - Normalized values.

`_combinatoril.txt` :
The file contains all sequences with processed intensity signals sorted from highest to lowest.
This is the format used for the original Seed-and-wobble Perl script.
The combinatorial file includes only probes with sequences provided in the fasta file.

`_masliner.png` $^*$ :
Figure describing the output of Masliner.

`_normalize.png` $^†$ :
Figure describing the array intensity layout according to the signal provided to the normalization and the normalized signal.

$^*$ - if Masliner was performed.

$^†$ - if normalization was performed.



### Example

```
python process.py -g1 './data/488nm_800_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr' -g2 './data/488nm_1000_80_1-KLF3LC_2-KLF3HC_3-SP1LC_4-SP1HC_5-IRF1_13.6.24_2-5.gpr' -f 'data/086902_D_Fasta_20220417.txt' -i 'F488 Median' -o './data/KLFLC'
```

### Notes


If Masliner is used before normalization, normalization is over the Adj values.
If disabling Masliner (by providing a single GPR file), normalization is over the selected intensity column from the original GPR (supplied by the `i` argument).

One can use the script without Masliner and without normalization by providing a single GPR and setting the radius to 0.
In that case, the script only merges the fasta file to the GPR and provides the corresponding "processed" and combinatorial file.

## Masliner

## Normalization

![Alt text](figures/norm_ver1.png)

## Reference

Berger, M. F., & Bulyk, M. L. (2009). Universal protein-binding microarrays for the comprehensive characterization of the DNA-binding specificities of transcription factors. Nature protocols, 4(3), 393-411r