# Analysis Pipeline for HIV Isoform Xploration (APHIX)
### Overview
APHIX is a pipeline for generating corrected HIV-1 spliced reads from end-to-end long-read sequences (such as from MrHAMER 2.0, Iso-Seq, or polyA Direct cDNA sequencing pipelines). The pipeline accepts FASTQ-format sequence files as input and outputs detailed information for each corrected isoform, splice site usage information and 
isoform prevelence information. 
### Features
The pipeline performs the following steps:
- Reads are mapped to reference genome using [minimap2] [minimap2]
- Mapped reads are converted to GFF using [spliced_bam2gff] [spliced_bam2gff] 
- Reads are clustered by exon/intron structure using [cluster_gff] [pinfish]
- Clusters are polished into splice junction-corrected reads using [polish_clusters][pinfish]
- Polished clusters are re-mapped and then converted back to GFF format
- Reads are preliminarily assigned an isoform identity using [gffcompare][gffcompare]
- Isoform identy is confirmed and splice site usage and isoform ratios are calculated using [HIV_Isoform_Checker][HIV_isoform_checker]
******************
# Getting Started
**_This pipeline is curretly only available for x86-64 processors and has only been tested on Ubuntu 20.04.5 LTS._**
### Requirements
The following software packages must be installed prior to running:

-  [miniconda3](https://conda.io/miniconda.html) - please refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Installation
After installing miniconda3, install the pipeline as follows:
```bash
git clone https://github.com/JessicaA2019/APHIX.git 
```
Change to directory:
```bash
cd APHIX
```
Create conda environment with all dependencies:
```bash
conda env create -f environment.yml
```
Activate environment:
```bash
conda activate APHIX
```
Copy precompiled modules provided by APHIX in conda environment:
```bash
cp precompiled_files/* $CONDA_PREFIX/bin
```
Deactivate environment:
```bash
conda deactivate
```
### Testing
To test if the installation was successful run:
```bash
conda activate APHIX
snakemake --cores 4 --configfile config.yml
```
 
******************
# Usage
In order for the pipeline to run, the following files/folders must be in the working directory:
| File | Location in APHIX directory|
|-------|-------------|
| Snakefile | ./Snakefile |
| config.yml | ./config.yml |
| whole auxillary_scripts folder | ./auxillary_scripts/* |
> Note: This pipeline was designed for HIV-1 references without any major insertions or deletions. If your sequences have larger deletions (ex: deltaENV or delta Matrix) or insertions (ex: multiple markers such as eGFP), the -l/--lengthFS parameter may need to be adjusted for HIV_Isoform_Checker. See the [GitHub][HIV_isoform_checker] for more details.

### Inputs
To run the pipeline the following input files are required:
| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. NL43) |
| Long-read sequences | A single concatenated FASTQ file. |

### Config File
Before starting the pipeline, confirm the following varibles are correct in the config.yml file:
| Input | Description |
|-------|-------------|
| sample_name | The name that the output folders/files will use. |
| input_fastq | The location of the long-read sequences fastq file.|
| reference_fasta | The location of the reference genome fasta file.|
| reference_gff_name | The name for reference gff to be created.|
| threads | The number of threads you would like the pipeline to use.|
| input_is_fastq | If the input_fastq file is a fastq, this should be True. If the input_fastq file is a fasta, this should be set to False.|
| cluster_size | The minimum size of a cluster during clustering and polishing.|

### Parameters in Snakefile
If the parameters for any step in the pipeline need to be adjusted, they can be easily found in the "Optional parameters" section. Please refer to the manual for each tool linked above for more details on the availble parameters.

### Running the pipeline
Once the config.yml is updated with the correct information and the APHIX environment is activated, simply run the following command:
```bash
snakemake --cores 4 --configfile config.yml
```
`--cores` specifies how many CPU cores will be used by the pipeline.

### Outputs
 The main output files created by the pipeline are:
| Output | Description |
|--------|-------------|
| Detailed Isoform Information | A csv of the detailed information on each isoform found. |
| Splice Site Usage Information | A csv of the counts and percent usage for each  donor site, acceptor site, and pairwise combination. |
| Isoform Prevelence Information | A csv of the counts and percentage of total spliced isoforms for each isoform type. |

After the a pipeline analysis has completed, these files can be found at `{sample_name}_isoform_analysis/8_corrected_isoforms/` e.g. `APHIX_test_isoform_analysis/8_corrected_isoforms/`.

**************************
## Help
For issues or bugs, please report them on the [issues page][issues]. 

## License
MIT - Copyright (c) 2024 Jessica Lauren Albert



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [minimap2]: <https://github.com/lh3/minimap2>
   [spliced_bam2gff]: <https://github.com/nanoporetech/spliced_bam2gff>
   [pinfish]: <https://github.com/nanoporetech/pinfish>
   [gffcompare]: <https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>
   [HIV_isoform_checker]: <https://github.com/JessicaA2019/HIV_Isoform_Checker>
   [issues]: <https://github.com/JessicaA2019/APHIX/issues>

