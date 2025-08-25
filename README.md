# PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data 

Precise profiling of epigenomes is essential for better understanding chromatin biology and gene regulation. Cleavage Under Targets & Tagmentation (CUT&Tag) is an easy and low-cost epigenomic profiling technique that can be performed on a low number of cells and at the single-cell level. With its growing adoption, CUT&Tag datasets spanning diverse biological systems are rapidly accumulating in the field. CUT&Tag assays use the hyperactive transposase Tn5 for DNA tagmentation. Tn5’s preference toward accessible chromatin alters CUT&Tag sequence read distributions in the genome and introduces open chromatin bias that can confound downstream analysis, an issue more substantial in sparse single-cell data. We show that open chromatin bias extensively exists in published CUT&Tag datasets, including those generated with recently optimized high-salt protocols. To address this challange, we present PATTY (Propensity Analyzer for Tn5 Transposase Yielded bias), a comprehensive computational method that corrects open chromatin bias in CUT&Tag data by leveraging accompanying ATAC-seq. By integrating transcriptomic and epigenomic data using machine learning and integrative modeling, we demonstrate that PATTY enables accurate and robust detection of occupancy sites for both active and repressive histone modifications, including H3K27ac, H3K27me3, and H3K9me3, with experimental validation. We further develop a single-cell CUT&Tag analysis framework built on PATTY and show improved cell clustering when using bias-corrected single-cell CUT&Tag data compared to using uncorrected data. Beyond CUT&Tag, PATTY sets a foundation for further development of bias correction methods for improving data analysis for all Tn5-based high-throughput assays.


## 0. Introduction to the PATTY Package

**PATTY** is a computational tool designed to correct open chromatin bias in CUT&Tag data at both **bulk** and **single-cell** levels. It leverages a pre-trained logistic regression model, built using CUT&Tag data in the K562 cell line, to correct bias for specific histone modifications.

- **Bulk mode:** PATTY applies the correction model to genome-wide 200bp tiling bins, and generates a bias-corrected score ranging from 0 to 1. Higher scores indicate higher confidence of true histone mark occupancy, while lower scores reflect likely false-positive or background signals due to open chromatin bias.

- **Single-cell mode:** PATTY performs bias correction at the individual cell level, producing a 200bp-bin by cell matrix of corrected signals. It then supports downstream cell clustering analysis using the bias-corrected data to improve biological interpretability and resolution.


- Changelog<br>
v1.0.0 First version of PATTY with both bulk and single-cell(sc) mode.

## 1. Installation
- Package requirements<br>
PATTY requires [Python](https://www.python.org) 3.6+ to run.<br>
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
- PATTY requires python3 packages [scipy](https://scipy.org) and [numpy](https://numpy.org) pre-installed. 
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
$ PATTY -m sc -c ${path}/testdata_scCUTTAGreads.bed.gz -a ${path}/testdata_ATACreads.bed.gz -f H3K27me3 -o testsc  
```

\# bulk mode 
```sh
$ SELMA -m bulk -c ${path}/testdata_bulkCUTTAGreads.bed.gz -a ${path}/testdata_ATACreads.bed.gz -f H3K27me3 -o testbulk 
```

## 4. Customize candidate peaks/regions.
PATTY provides an option (-p) to take user-supplied customized peak files as the target regions for the PATTY analysis. The required peak file should be in BED format (plain text), have >=4 columns (chrom, start, end, name), and the total length should be >= 200kbps (1k 200bp bins). By default, PATTY uses the peaks detected from the same dataset (e.g., fragments.bed file) using SICER for the peak calling (with any method) to ensure sufficient cleavages/signal on the peak regions. Below is the example with an external/customized peak file: <br />
The test files (testdata_reads.bed.gz and testpeak.bed) in the following cmd lines can be downloaded via the link in section x

\# sc mode 
```sh
$ PATTY -m sc -c ${path}/testdata_CUTTAGreads.bed.gz -a ${path}/testdata_ATACreads.bed.gz -p ${path}/testpeak.bed  -f H3K27me3 -o testsc  
```

\# bulk mode 
```sh
$ SELMA -m bulk -c ${path}/testdata_CUTTAGreads.bed.gz -a ${path}/testdata_ATACreads.bed.gz -p ${path}/testpeak.bed  -f H3K27me3 -o testbulk 
```

## 5. Pre-processing Steps for Generating the Input Fragments File

PATTY takes aligned fragment files in **BED format** as input. Users may apply any preferred pre-processing pipeline to generate these files. We recommend retaining only **high-quality reads** with **MAPQ > 30** to ensure accurate bias correction.

### Default Input Format

The expected BED format varies depending on data type:

#### • Bulk CUT&Tag – Paired-End
Paired-end data should be pre-processed into fragment-level BED format:

- `chr1`, `10500`, `10646`: Chromosome, start, and end of the fragment
- `.`: Placeholder for name
- `60`: Mapping quality (optional)
- `+` or `-`: Strand (optional)

> Ensure duplicates are removed and only **unique fragments** are retained.

#### • Bulk CUT&Tag – Single-End
For single-end reads, convert each read into a **fixed-length fragment** (e.g., 146 bp) extending from the 5′ end:



## 6. Install and use published single-cell clustering methods based on PATTY bias correction. 
PATTY sc mode implements several cell clustering methods in the single-cell clustering analysis in addition to the default Kmeans analysis. To activate these methods (name, version and link listed below), users need to install the related package, and specify the method by the --clusterMethod parameter. If a method is declared by the --clusterMethod parameter but is not installed, SELMA will skip the single-cell clustering analysis.
- Seurat (required packages: [ArchR v1.0.1](https://satijalab.org/seurat/), [tabix](http://www.htslib.org/doc/tabix.html), and [bgzip](http://www.htslib.org/doc/bgzip.html))
- scran (required packages: [ArchR v1.0.1](https://satijalab.org/seurat/), [tabix](http://www.htslib.org/doc/tabix.html), and [bgzip](http://www.htslib.org/doc/bgzip.html))
- APEC (required package: [APEC v1.2.2](https://github.com/QuKunLab/APEC))

PATTY also provides UMAP/t-SNE visualization for the single-cell clustering analysis. Users can activate this function by the --UMAP parameter. For the PCAkm method, the [umap](https://cran.r-project.org/web/packages/umap/index.html) package in R is required. 

## 7. Output files
1. `NAME_summaryReports.pdf` is the summary pdf file which contains information on:
     - Input file and parameter description
     - basic QC of the data
     - Summary of the SELMA bias estimation/correction results

    \#Note: This pdf file is only generated if pdflatex is pre-installed. A NAME_summaryReports.txt file is generated as well. A .tex file will also be generated in case users want to make the pdf document later.

2. `NAME_peaks.bed` is the peaks detected from the fragment files (using SICER). Each peak was split into 200bp bins. 

3. `NAME_PATTYscore.bdg` (bulk mode only) is the PATTYscore for each candidate 200bp-bin, in bedGraph format. The scores are transformed into a 200bp-resolution track for the input regions. For genomic regions not covered by the input regions, PATTY will assign 0 as the score. 

4. `NAME_binXcell.txt.gz` (sc mode only) is the bin-by-cell PATTY score matrix generated from the single-cell analysis. Cells are filtered by the total reads count per cell (default >=10,000 reads, can be changed through the parameters). 

5. `NAME_scClustering.txt.gz` (sc mode only) is the cell clustering result using the PATTY bias corrected matrix.


## 8. Testing data and example of output files
We provided the test data for users to test SELMA. The sc/bulk output can also be generated with the cmd lines in Section 3/4 using the testing data as input. Click the file names to download (copy the backupLink for cmdline download). 
- testing data: [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_reads.bed.gz)
- testing peak file(optional for -p): [`Dropbox`](https://www.dropbox.com/s/a4r3gzux7v72rr9/testpeak.bed)
- testing cellnames (optional for --cellnames in sc mode): [`Dropbox`](https://www.dropbox.com/s/9eb60r9xx7gbh13/testsc_cellnames.txt)
- output for PATTY **bulk** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/x8f29ao73t5ka8a/AADPjRgtgmW0DXJTiPMYWIS-a?dl=0)
- output for PATTY **sc** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/jl7a28w9a984tpz/AACoGYzLRnBwZLgmbGlu6bwSa?dl=0) 
- The PATTY with testing data (e.g., using sc mode) will be finished within 30 minutes.


## 9. Other parameters in the PATTY pipeline
You can also set the following parameters for more accurate bias estimation and correction:
- -\-cellnames=CELLNAMES  
[sc optional] Single column file for name list of used individual cells, each line contains the name of an individual cell. This parameter is only used for sc mode. This parameter is not used very commonly. 
- -\-readCutoff=READCUTOFF  
[sc optional] Reads number cutoff for high-quality cells. Cells with < 10000(default) reads will be discarded in the analysis. Users can change this parameter for samples with low sequencing depth to include more cells in the analysis. Setting a lower number for this parameter will possibly decrease the accuracy of clustering results due to the low-quality cells. 
- -\-peakMinReads=PEAKMINREADS  
[sc optional] Peaks with < 10(default) cleavages covered (across all high-quality cells) will be discarded in the analysis.
- -\-peakMaxReads=PEAKMAXREADS  
[sc optional] Peaks with > X cleavages covered (across all high-quality cells) will be discarded in the analysis. Set 0 to close this function (default)
- -\-clusterMethod=CLUSTERMETHOD  
[sc optional] Method used for single-cell clustering analysis. The default is Kmeans(PCA dim reduction + K-means clustering). Optional choices (Seurat, scran, and APEC) require related packages installed (described in section x)
- -\-clusterNum=CLUSTERNUM  
[sc optional] Number of clusters specified for K-means clustering and only used for the PCAkm (setting by --clusterMethod) method. The default is 7. 
- -\-topDim=TOPDIM  
[sc optional] Number of dimensions (with highest Variance) used for clustering. Only used for PCAkm(PC) and ArchR (Latent variable). This number is suggested to be >=30 (deafult=60)
- -\-UMAP  
[sc optional] Turn on this parameter to generate a UMAP plot for the clustering results
- -\-overwrite  
[optional] Force overwrite; setting this parameter will remove the existing result! PATTY will terminate if there is a folder with the same name as -o in the working directory. Set this parameter to force PATTY to run. 
- -\-keeptmp  
[optional] Whether or not to keep the intermediate results (tmpResults/)

# Reproduce cell clustering results using the PATTY package
Users can reproduce the clustering results for the nano-CT data in the manuscript (Figure 6, H3K27me3 nano-CT data in mouse brain,  K-means clustering) by running PATTY with the following cmd line:
```sh
$ PATTY -m sc -i ${path}/testdata_reads.bed.gz -f H3K27me3 -o nanoCT_H3K27me3 --UMAP --overwrite --keeptmp --cellnames ${path}/testsc_cellnames.txt
```
The test files (testdata_reads.bed.gz and testsc_cellnames.txt) in the cmd line can be downloaded via the link in section x.

# Supplementary data
- PATTY pre-trained models for [H3K27me3](https://www.dropbox.com/s/ncemdhp0cee3cic/DNase_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/DNase_SELMAbias_10mer.txt.gz)) and [H3K27ac](https://www.dropbox.com/s/x5iiy27ef80fl19/ATAC_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/ATAC_SELMAbias_10mer.txt.gz)). Both models were trained from CUT&Tag data in K562 cell line. Note that these pre-trained models are already built-in and used in the PATTY package.  

