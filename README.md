# PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data 

Precise profiling of epigenomes is essential for better understanding chromatin biology and gene regulation. Cleavage Under Targets & Tagmentation (CUT&Tag) is an easy and low-cost epigenomic profiling technique that can be performed on a low number of cells and at the single-cell level. With its growing adoption, CUT&Tag datasets spanning diverse biological systems are rapidly accumulating in the field. CUT&Tag assays use the hyperactive transposase Tn5 for DNA tagmentation. Tn5’s preference toward accessible chromatin alters CUT&Tag sequence read distributions in the genome and introduces open chromatin bias that can confound downstream analysis, an issue more substantial in sparse single-cell data. We show that open chromatin bias extensively exists in published CUT&Tag datasets, including those generated with recently optimized high-salt protocols. To address this challange, we present PATTY (Propensity Analyzer for Tn5 Transposase Yielded bias), a comprehensive computational method that corrects open chromatin bias in CUT&Tag data by leveraging accompanying ATAC-seq. By integrating transcriptomic and epigenomic data using machine learning and integrative modeling, we demonstrate that PATTY enables accurate and robust detection of occupancy sites for both active and repressive histone modifications, including H3K27ac, H3K27me3, and H3K9me3, with experimental validation. We further develop a single-cell CUT&Tag analysis framework built on PATTY and show improved cell clustering when using bias-corrected single-cell CUT&Tag data compared to using uncorrected data. Beyond CUT&Tag, PATTY sets a foundation for further development of bias correction methods for improving data analysis for all Tn5-based high-throughput assays.


## 0. Introduction to the PATTY Package

**PATTY** is a computational tool designed to correct open chromatin bias in CUT&Tag data at both **bulk** and **single-cell** levels. It leverages a pre-trained logistic regression model, built using CUT&Tag data in the K562 cell line, to correct bias for specific histone modifications.

- **Bulk mode:** PATTY applies the correction model to genome-wide 200bp tiling bins, and generates a bias-corrected score on each candidate bin. 

- **Single-cell mode:** PATTY performs bias correction at the individual cell level, producing a 200bp-bin by cell matrix of corrected signals. It then supports downstream cell clustering analysis using the bias-corrected data to improve biological interpretability and resolution.


- Changelog<br>
v1.0.0 First version of PATTY with both bulk and single-cell(sc) mode.

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
- Bedtools (Quinlan et al., Bioinformatics, 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) are recommended if users want to generate bigwig tracks for bias-corrected data. 
- The installation should be finished in about one minute.


## 2. Run PATTY (usage)
#### Essential parameters
To run PATTY with the default parameters, you can set the following parameters:
-   -m MODE, --mode=MODE
Mode of PATTY, choose from sc(single-cell) or bulk
-   -c CUTTAG, --cuttag=CUTTAG
Input fragments file in (paired/single end) bed format for CUT&Tag data, with .bed extension. For sc mode, the 4th(name) column of the bed file represents the name/barcode of the corresponding individual cell.
-   -a ATAC, --atac=ATAC
Input fragments file in bed format for ATAC-seq data, with .bed extension. The ATAC-seq fragments were used as bulk data for both sc and bulk modes. 
-   -f FACTOR, --factor=FACTOR
Factor type of the CUT&Tag data. Currently PATTY support H3K27me3 (default), H3K27ac, and H3K9me3
-   -o OUTNAME, --outname=OUTNAME
Name of output results

Example of running SELMA with default parameters (test data downloadable in :

\# sc mode 
```sh
$ PATTY -m sc -c ${path}/testdata_scCUTTAGreads.bed.gz -a ${path}/testdata_scATACreads.bed.gz -f H3K27me3 -o testsc  
```

\# bulk mode 
```sh
$ SELMA -m bulk -c ${path}/testdata_bulkCUTTAGreads.bed.gz -a ${path}/testdata_bulkATACreads.bed.gz -f H3K27me3 -o testbulk 
```

## 4. Pre-processing Steps for Generating the Input Fragments File

PATTY takes aligned fragment files in **BED format** as input. Users may apply any preferred pre-processing pipeline to generate these files. We recommend retaining only **high-quality reads** with **MAPQ > 30** to ensure accurate bias correction. Note that PATTY takes original fragments bed files as input (e.g., transformed directly from aligned BAM files, or 10x cell ranger outputed fragments.tsv file for sc data). Please don't do any customized extension or shifting. 

### Default Input Format

The expected BED format varies depending on data type:

#### • Bulk CUT&Tag (Paired-End)
```
chr1    10500   10646
chr2    20840   21000
```

#### • Bulk CUT&Tag (Single-End)
```
chr1    10500   10646   .   .   +
chr2    20840   20986   .   .   -
```
> The 4-5th column represents an optional placeholder.

#### • Single-Cell CUT&Tag
```
chr1    10500   10646   CellA
chr2    20840   21000   CellB
```
> The 4th column must contain the cell barcode or cell name (like AATAACTACGCC-1).


## 5. Install and use published single-cell clustering methods based on PATTY bias correction. 
PATTY sc mode implements several cell clustering methods in the single-cell clustering analysis, in addition to the default K-means analysis. To activate these methods, users need to install the related package and specify the method by the --clusterMethod parameter. If a method is declared by the --clusterMethod parameter but is not installed, SELMA will skip the single-cell clustering analysis.

PATTY also provides UMAP visualization for the single-cell clustering analysis. Users can activate this function by using the --UMAP parameter. For this option, the [umap](https://cran.r-project.org/web/packages/umap/index.html) package in R is required. 

## 6. Output Files

### Bulk Mode Outputs

1. `NAME_PATTYscore.bdg`  
   A 200bp-resolution genome-wide track in **bedGraph** format containing the PATTY scores for each candidate bin.  
   - Scores range from 0 to 1. Higher scores indicate higher confidence of true histone mark occupancy, while lower scores reflect likely false-positive or background signals due to open chromatin bias.

### Single-Cell Mode Outputs

1. `NAME_binXcell.txt.gz`  
   A **bin-by-cell PATTY score matrix** generated from single-cell CUT&Tag analysis.  
   - Rows: 200bp bins  
   - Columns: individual cells  

2. `NAME_scClustering.txt.gz`  
   The cell clustering result is based on the PATTY bias-corrected matrix.  
   - Format: tab-delimited, with each cell's cluster label


## 9. Testing data and example of output files
We provided the test data for users to test PATTY. The sc/bulk output can also be generated with the command lines in Section 2 using the testing data as input. Click the file names to download. 
- testing data for **bulk** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_bulkCUTTAGreads.bed.gz)
   - ATAC [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_bulkATACreads.bed.gz)
- testing data for **sc** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_scCUTTAGreads.bed.gz)
   - ATAC [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_scATACreads.bed.gz)
- output for PATTY **bulk** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/x8f29ao73t5ka8a/AADPjRgtgmW0DXJTiPMYWIS-a?dl=0)
- output for PATTY **sc** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/jl7a28w9a984tpz/AACoGYzLRnBwZLgmbGlu6bwSa?dl=0) 
- The PATTY with testing data (e.g., using sc mode) will be finished within 60 minutes.


## 9. Other parameters in the PATTY pipeline
You can also set the following parameters for more accurate bias estimation and correction:
- -\-binMinReads=BINMINREADS  
[optional] Bins with < 5(default) reads covered will be discarded in the analysis. For sc mode, bins with a total of < 5 (default) reads across all high-quality cells will be discarded. 
- -\-readCutoff=READCUTOFF  
[sc optional] Reads number cutoff for high-quality cells. Cells with < 10000(default) reads will be discarded in the analysis. Users can change this parameter for samples with low sequencing depth to include more cells in the analysis. Setting a lower number for this parameter may decrease the accuracy of clustering results due to the low-quality cells. 
- -\-binMaxReads=BINMAXREADS  
[sc optional] Bins with > X cleavages covered (across all high-quality cells) will be discarded in the analysis. Set 0 to close this function (default)
- -\-clusterMethod=CLUSTERMETHOD  
[sc optional] Method used for single-cell clustering analysis. The default is K-means (PCA dim reduction + K-means clustering). Optional choices (Seurat and scran) require related packages installed (described in section x)
- -\-clusterNum=CLUSTERNUM  
[sc optional] Number of clusters specified for K-means clustering and only used for the PCAkm (setting by --clusterMethod) method. The default is 7. 
- -\-UMAP  
[sc optional] Turn on this parameter to generate a UMAP plot for the clustering results
- -\-overwrite  
[optional] Force overwrite; setting this parameter will remove the existing result! PATTY will terminate if there is a folder with the same name as -o in the working directory. Set this parameter to force PATTY to run. 
- -\-keeptmp  
[optional] Whether or not to keep the intermediate results (tmpResults/)

# Reproduce cell clustering results using the PATTY package
Users can reproduce the bulk correction data in the manuscript (Figure 4, H3K27me3 CUT&Tag data in K562) by running PATTY with the following command line:
```sh
$ PATTY -m bulk -c ${path}/testdata_bulkCUTTAGreads.bed.gz -a ${path}/testdata_bulkATACreads.bed.gz -f H3K27me3 -o K562_H3K27me3 --UMAP --overwrite --keeptmp 
```
Users can reproduce the clustering results for the nano-CT data in the manuscript (Figure 7, H3K27me3 nano-CT data in mouse brain,  K-means clustering) by running PATTY with the following command line:
```sh
$ PATTY -m sc -c ${path}/testdata_scCUTTAGreads.bed.gz -a ${path}/testdata_scATACreads.bed.gz -f H3K27me3 -o nanoCT_H3K27me3 --UMAP --overwrite --keeptmp 
```
The test data in the command line can be downloaded via the link in section 9.


