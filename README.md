# spacerAQ_vh
Custom for analyzing data in Hoikkala et al. 2020. Currently for review purposes only. The pipeline has been tested with OSX 10.14.6.

# Instructions
The analysis can be run using the pipeline spacerpipeline_hoikkala_etal_2020.sh. However, several steps need to be taken before the run:

## 1. Install dependencies
The list of necessary dependencies are in the file environment.yml. The easiest way to install them is to create a new Anaconda environment with the command `conda env create -f environment.yml` while in the root folder. This will create a new conda environment named spacerAQ_vh_2020

## 2. Place trimmomatic-0.36.jar in root folder
Or modify the script to indicate where your trimmomatic file is located. It is also possible to skip trimming - this will not dramatically alter the results (but requires you to modify file paths in spacerpipeline_hoikkala_etal_2020.sh.

