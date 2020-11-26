# spacerAQ_vh
Pipeline for analyzing data in Hoikkala et al. 2020. Currently for review purposes only.

Typical install time: 10 min.\n
Typical run time of the pipeline: 2-5 minutes without simulations, one hour with simulations (simulations currently disabled; to enable, uncomment line 264 from spacerpipeline_absolute_hoikkala_etal_2020.sh or from spacerpipeline_unique_hoikkala_etal_2020.sh).
Required operating system: Unix (tested with MacOS 10.15.7 with zsh).
Output: several files and folders regarding extracted spacers and their properties. Running the pipeline with default settings will reproduce the results presented in the manuscript. R scripts can be used to recreate figures of the manuscript after data has been generated. Separate pipeline scripts for running analyses on either absolute or unique spacers (spacerpipeline_absolute_hoikkala_etal_2020.sh or spacerpipeline_unique_hoikkala_etal_2020.sh)

# Instructions
The analysis is run using the shell pipeline spacerpipeline_absolute_hoikkala_etal_2020.sh or spacerpipeline_unique_hoikkala_etal_2020.sh.
A little setup is required after downloading the files on your computer (remember to maintain folder structure)

## 1. Install dependencies
The list of necessary dependencies are in the file environment.yml. The easiest way to install them is to create a new [Anaconda](https://www.anaconda.com/distribution/) environment with the command `conda env create -f environment.yml` while in the root folder. This will create a new conda environment named spacerAQ_vh_2020. Activate this environment to continue.

## 2. Download raw data files
Link for the raw data is in the materials and methods file of the manuscript. Copy these fastq files to the folder pipeline/raw_reads.

## 3. Run the pipeline
Run the pipeline `sh spacerpipeline_absolute_hoikkala_etal_2020.sh` or `sh spacerpipeline_unique_hoikkala_etal_2020.sh`

# R
Once the pipeline has been run, the provided R files can be used to illustrate the results.

For any questions, create a new github issue or email me at ville.hoikkala@jyu.fi
