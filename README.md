# spacerAQ_vh
Pipeline for analyzing data in Hoikkala et al. 2020 (submitted to Nature). Currently for review purposes only.

Typical install time: 10 min. Typical run time of the pipeline: 2-5 minutes without simulations, one hour with simulations (simulations currently disabled; to enable, uncomment line 264 from spacerpipeline_hoikkala_etal_2020.sh). Required operating system: Unix (tested with MacOS 10.14.6 with shell bash). Output: several files and folders regarding extracted spacers and their properties. Running the pipeline with default settings will reproduce the results presented in the manuscript. R scripts can be used to recreate figures of the manuscript after data has been generated.

# Instructions
The analysis is run using the shell pipeline spacerpipeline_hoikkala_etal_2020.sh.
Several steps need to be taken to setup the run after downloading the files on your computer (remember to maintain folder structure)

## 1. Install dependencies
The list of necessary dependencies are in the file environment.yml. The easiest way to install them is to create a new [Anaconda](https://www.anaconda.com/distribution/) environment with the command `conda env create -f environment.yml` while in the root folder. This will create a new conda environment named spacerAQ_vh_2020. Activate this environment to continue.

## 2. Locate your trimmomatic-0.36.jar and copy it in the root folder (same folder where environment.yml is located)
... or modify the script to indicate where your trimmomatic file is located. It is also possible to skip trimming - this will not dramatically alter the results (but requires you to modify file paths in spacerpipeline_hoikkala_etal_2020.sh).

## 3. Download raw data files
Link for the raw data is in the materials and methods file. Copy these fastq files to the folder pipeline/raw_reads.

## 4. Run the pipeline
Use command `sh spacerpipeline_hoikkala_etal_2020.sh`

# R
Once the pipeline has been run, the provided R files can be used to illustrate the results
