# PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data 

Precise profiling of epigenomes is essential for better understanding chromatin biology and gene regulation. Cleavage Under Targets & Tagmentation (CUT&Tag) is an easy and low-cost epigenomic profiling technique that can be performed on a low number of cells and at the single-cell level. With its growing adoption, CUT&Tag datasets spanning diverse biological systems are rapidly accumulating in the field. CUT&Tag assays use the hyperactive transposase Tn5 for DNA tagmentation. Tn5’s preference toward accessible chromatin alters CUT&Tag sequence read distributions in the genome and introduces open chromatin bias that can confound downstream analysis, an issue more substantial in sparse single-cell data. We show that open chromatin bias extensively exists in published CUT&Tag datasets, including those generated with recently optimized high-salt protocols. To address this challange, we present PATTY (Propensity Analyzer for Tn5 Transposase Yielded bias), a comprehensive computational method that corrects open chromatin bias in CUT&Tag data by leveraging accompanying ATAC-seq. By integrating transcriptomic and epigenomic data using machine learning and integrative modeling, we demonstrate that PATTY enables accurate and robust detection of occupancy sites for both active and repressive histone modifications, including H3K27ac, H3K27me3, and H3K9me3, with experimental validation. We further develop a single-cell CUT&Tag analysis framework built on PATTY and show improved cell clustering when using bias-corrected single-cell CUT&Tag data compared to using uncorrected data. Beyond CUT&Tag, PATTY sets a foundation for further development of bias correction methods for improving data analysis for all Tn5-based high-throughput assays.


## 0. Introduction to the PATTY Package

**PATTY** is a computational tool designed to correct open chromatin bias in CUT&Tag data at both **bulk** and **single-cell** levels. Current version of PATTY support open chromatin bias correction for H3K27me3, H3K27ac, and H3K9me3. It leverages a pre-trained logistic regression model, built using CUT&Tag data in the K562 cell line, to correct bias for specific histone modifications.

- **Bulk mode:** PATTY applies the correction model to genome-wide 200bp tiling bins, and generates a bias-corrected score on each candidate bin. 

- **Single-cell mode:** PATTY performs bias correction at the individual cell level, producing a 200bp-bin by cell matrix of corrected signals. It then supports downstream cell clustering analysis using the bias-corrected data to improve biological interpretability and resolution.


- Changelog<br>
v1.0.0 PATTY for biorxiv manuscript and initial submission 

## 1. Installation
- Package requirements<br>
PATTY requires [Python](https://www.python.org) 3.6+ and [Rscript](https://www.r-project.org) v3+ to run.<br>
PATTY requires Python packages [scipy](https://scipy.org) and [numpy](https://numpy.org) pre-installed.

\# for root user
```sh
$ cd PATTY
$ sudo python setup.py install  
```
\# if you are not the root user, you can install PATTY at a specific location where you have write permission
```sh
$ python setup.py install --prefix /home/PATTY  # Here you can replace “/home/PATTY” with any location 
$ export PATH=/home/PATTY/bin:$PATH    # setup PATH for the software
$ export PYTHONPATH=/home/PATTY/lib/python3.6/site-packages:$PYTHONPATH    # setup PYTHONPATH for module import
```
\# To check the PATTY package, just type:
```sh
$ PATTY --help  # If you see the help manual, you have successfully installed PATTY
```

\# NOTE: 
- To install PATTY on MacOS, the users need to download and install Command Line Tools beforehand
- Bedtools (Quinlan et al., Bioinformatics, 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) will be installed automatically if not installed. 

## 2. Run PATTY (usage)
#### Essential parameters
To run PATTY with the default parameters, you can set the following parameters:
-   -m MODE, --mode=MODE
Mode of PATTY, choose from bulk or sc(single-cell)
-   -c CUTTAG, --cuttag=CUTTAG
CUTTAG Input fragments file in (paired/single end) bed format for CUT&Tag data, with .bed extension (or .bed.gz for comparessed file). For sc mode, the 4th(name) column of the bed file represents the name/barcode of the corresponding individual cell
-   -a ATAC, --atac=ATAC
ATAC Input fragments file in bed format for ATAC-seq data, with .bed extension(or .bed.gz for comparessed file). The ATAC-seq fragments were used as bulk data for both sc and bulk modes (only chrm,start,end 3columns are required)
-   -f FACTOR, --factor=FACTOR
FACTOR Factor type of the CUT&Tag data. Currently PATTY support H3K27me3 (default), H3K27ac, and H3K9me3
-   -g GENOME, --genome=GENOME
genome version of the input data, choose from hg38 (default) and mm10
-   -o OUTNAME, --outname=OUTNAME
Name of output results

Example of running PATTY with default parameters (test data downloadable in :

\# bulk mode 
```sh
$ PATTY -m bulk -c ${path}/testdata_bulk_H3K27me3_reads.bed.gz -a ${path}/testdata_bulk_ATAC_reads.bed.gz -f H3K27me3 -o testbulk 
```

\# sc mode 
```sh
$ PATTY -m sc -c ${path}/testdata_sc_H3K27me3_reads.bed.gz -a ${path}/testdata_sc_ATAC_reads.bed.gz -f H3K27me3 -o testsc  
```


## 4. Pre-processing Steps for Generating the Input Fragments File

PATTY takes aligned fragment files in **BED format** as input(or .bed.gz for gzip comparessed file). Users may apply any preferred pre-processing pipeline to generate these files. We recommend retaining only **high-quality reads** with **MAPQ > 30** to ensure accurate bias correction. Note that PATTY takes original fragments bed files as input (e.g., transformed directly from aligned BAM files, or 10x cell ranger outputed fragments.tsv file for sc data). Please don't do any customized extension or shifting. 

### Default Input Format

The expected BED format varies depending on data type:

#### • Bulk CUT&Tag (Single-End)
```
chr1    10500   10646   .   .   +
chr2    20840   20986   .   .   -
```
> The 4-5th column represents an optional placeholder.

#### • Bulk CUT&Tag (Paired-End)
```
chr1    10500   10646
chr2    20840   21000
```

#### • Single-Cell CUT&Tag
```
chr1    10500   10646   CellA
chr2    20840   21000   CellB
```
> The 4th column must contain the cell barcode or cell name (like AATAACTACGCC-1).

## 5. Output Files

### Bulk Mode Outputs

1. `NAME_PATTYscore.bw`  
   A 200bp-resolution genome-wide track in **bigWig** format containing the PATTY scores for each candidate bin.  
   - Scores range from 0 to 1. Higher scores indicate higher confidence of true histone mark occupancy, while lower scores reflect likely false-positive or background signals due to open chromatin bias.

### Single-Cell Mode Outputs

1. `NAME_binXcell.txt.gz`  
   A **bin-by-cell PATTY score matrix** generated from single-cell CUT&Tag analysis.  
   - Rows: 200bp bins  
   - Columns: individual cells
   - Values: Similar PATTY score like in Bulk mode but for each individual cell  

2. `NAME_scClustering.txt.gz`  
   The cell clustering result is based on the PATTY bias-corrected matrix.  
   - Format: tab-delimited, with each cell's cluster label


## 6. Testing data and example of output files
We provided the test data for users to test PATTY. The sc/bulk output can also be generated with the command lines in Section 2 using the testing data as input. Click the file names to download. 
- testing data for **bulk** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/scl/fi/820c3ryhj7ffbbkysgd5k/testdata_bulk_H3K27me3_reads.bed.gz)
   - ATAC [`Dropbox`](https://www.dropbox.com/scl/fi/n90qbq6xf1hir31legnlc/testdata_bulk_ATAC_reads.bed.gz)
- testing data for **sc** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/scl/fi/t9hv9okgubvmevbafysh5/testdata_sc_H3K27me3_reads.bed.gz)
   - ATAC [`Dropbox`](https://www.dropbox.com/scl/fi/tlb9tzn32ykwzlh19znch/testdata_sc_ATAC_reads.bed.gz)
- output for PATTY **bulk** mode using bulk-testing data input: [`Dropbox`](https://www.dropbox.com/scl/fi/597k0encidxcxr02m5db7/testdata_bulk_H3K27me3_PATTYcorrect.bw)
- output for PATTY **sc** mode using sc-testing data input: [`Dropbox`](https://www.dropbox.com/scl/fi/8t5cjaani27tgs82lfpcr/testdata_sc_H3K27me3_clustering.txt.gz) 
- The PATTY with testing data (e.g., using sc mode) will be finished within 60 minutes.


## 7. Other parameters in the PATTY pipeline
You can also set the following parameters for more accurate bias estimation and correction:
- -\-binMinReads=BINMINREADS  
[optional] Bins with < 5(default) reads covered will be discarded in the analysis. For sc mode, bins with a total of < 5 (default) reads across all high-quality cells will be discarded. set 0 to turn off this parameter. 
- -\-binList=BINLIST  
[optional] Bed file for inputting candidate bins/peaks for the analysis. When inputted, the correction will be done in only these bins for bulk mode. The correction/clustering will be only done these bins for sc mode. This parameter is designed for customized high-reads bin (bulk mode) or high-var bin (sc mode). The inputed peaks/bins will be transformed/splited to 200bp bins for as the input. 
- -\-cellnames=CELLNAMES  
[optional] Single column plain text file for name list of used individual cells, each line contain the name of the individual cell. This parameter is only used for sc mode. 
- -\-readCutoff=READCUTOFF  
[sc optional] Reads number cutoff for high-quality cells. Cells with < 10000(default) reads will be discarded in the analysis. Users can change this parameter for samples with low sequencing depth to include more cells in the analysis. Setting a lower number for this parameter may decrease the accuracy of clustering results due to the low-quality cells. 
- -\-clusterMethod=CLUSTERMETHOD  
[sc optional] Method used for single-cell clustering analysis. The default is K-means (PCA dim reduction + K-means clustering). Optional choices (Seurat and scran) require related packages installed (described in section x)
- -\-clusterNum=CLUSTERNUM  
[sc optional] Number of clusters specified for K-means clustering and only used for the PCAkm (setting by --clusterMethod) method. The default is 7. 
- -\-UMAP  
[sc optional] Turn on this parameter to generate a UMAP plot for the clustering results.
- -\-overwrite  
[optional] Force overwrite; setting this parameter will remove the existing result! PATTY will terminate if there is a folder with the same name as -o in the working directory. Set this parameter to force PATTY to run. 
- -\-keeptmp  
[optional] Whether or not to keep the intermediate results (tmpResults/)

## 8. Reproduce cell clustering results using the PATTY package
Users can reproduce the bulk correction data in the manuscript (Figure 4A, H3K27me3 CUT&Tag data in K562) by running PATTY with the following command line:
```sh
$ PATTY -m bulk -c ${path}/testdata_bulk_H3K27me3_reads.bed.gz -a ${path}/testdata_bulk_ATAC_reads.bed.gz -f H3K27me3 -o testdata_bulk_H3K27me3 
```
Users can reproduce the clustering results for the sc correction data in the manuscript (Figure 7B, H3K27me3 nano-CT data in mouse brain,  K-means clustering) by running PATTY with the following command line:
```sh
$ PATTY -m sc -c ${path}/testdata_sc_H3K27me3_reads.bed.gz -a ${path}/testdata_sc_ATAC_reads.bed.gz -f H3K27me3 -o testdata_sc_H3K27me3 --UMAP --overwrite --keeptmp 
```
The example input data in the command line and output results can be downloaded via the link in section 6.


