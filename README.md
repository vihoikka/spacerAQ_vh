# spacerAQ_vh
Pipeline for analyzing data in Hoikkala et al. 2020: "Shared CRISPR-Cas spacer acquisition between subtypes VI-B and II-C" (submitted to Nature). Currently for review purposes only. The pipeline has been tested with OSX 10.14.6 using Shell bash.

# Instructions
The analysis can be run using the pipeline spacerpipeline_hoikkala_etal_2020.sh.
Several steps need to be taken to setup the run:

## 1. Install dependencies
The list of necessary dependencies are in the file environment.yml. The easiest way to install them is to create a new [Anaconda](https://www.anaconda.com/distribution/) environment with the command `conda env create -f environment.yml` while in the root folder. This will create a new conda environment named spacerAQ_vh_2020. Activate this environment to continue.

## 2. Place your trimmomatic-0.36.jar in root folder (if the current .jar does not work)
... or modify the script to indicate where your trimmomatic file is located. It is also possible to skip trimming - this will not dramatically alter the results (but requires you to modify file paths in spacerpipeline_hoikkala_etal_2020.sh).

## 3. Download raw data files
Files are in https://bit.ly/2tVWGhW. Copy these files to the folder pipeline/raw_reads.

## 4. Run the pipeline
Use command `sh spacerpipeline_hoikkala_etal_2020.sh`

# R
Once the pipeline has been run, the provided R files can be used to illustrate the data
