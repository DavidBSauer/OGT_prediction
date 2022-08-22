# OGT_prediction
Scripts for calculating features from genomic sequences, multiple linear regression of those features to the originating species' Optimal Growth Temperatures (OGT), and prediction of a species' OGT using those linear regression models.

See: Sauer & Wang. Predicting the optimal growth temperatures of prokaryotes using only genome derived features. Bioinformatics (2019) https://doi.org/10.1093/bioinformatics/btz059

## Installation and Requirements
This has been developed and tested on Ubuntu 18.04 LTS. The scripts should work on any system supporting Python 3, so long as the external programs are installed properly.

1. Download these scripts. This is easiest using git to clone the repository.
```
git clone https://github.com/DavidBSauer/OGT_prediction
```

2. Install the requirements.
These scripts depend upon the programs: Python3, tRNAscan-SE, Bedtools, Barrnap, and Prodigal. These have their own dependencies also.
The following python packages also need to be installed: numpy, scipy, matplotlib, biopython, bcbio-gff, tqdm, sklearn, matplotlib, and matplotlib-venn.

To install everything in Ubuntu (or other system that use the apt package manager), go into the downloaded directory and use the pre-made bash script.
```
cd OGT_prediction
./Ubuntu_setup.bash
```
If you're on another OS, install all the python packages. Then install bedtools, tRNAscan-SE, barrnap, and prodigal; and all their dependencies. Create a file called external_tools.txt listing bedtools, tRNAscan-SE, barrnap, and prodigal tab separated from the absolute path for each executable. Copy this external_tools.txt file into the feature_calculation and prediction directories.

### OGT Data
To calculate regressions, a species_trait file (tab delimited) is required. The species_OGT file used for calculating regression in Sauer & Wang (2019) was originally published in Sauer et al. (2015). The species_OGT file is named "all_merged_12_10_2012.txt" and found within http://www.med.nyu.edu/skirball-lab/dwanglab/files/thermostability_v1.1.tar.gz. If new regressions are calcuated using that species-OGT file, please cite: Sauer et al. Rapid Bioinformatic Identification of Thermostabilizing Mutations. Biophysical Journal (2015) https://doi.org/10.1016/j.bpj.2015.07.026.

## Demonstration
This uses the previously computed regression models to predict OGTs for a few species. Start within the prediction directory.

0. Move into the prediction directory.
```
cd prediction
```

1. Download genomes species IN the provided list.
```
python3 genome_retriever.py ../data/prediction_demo/species.txt IN 40
```

2. Download taxonomic classification for species.
```
python3 clade_retriever.py ../data/prediction_demo/species.txt your_email_address@awesome.com
```

3. Run the prediction script to predict the OGT of each species.
```
python3 prediction_pipeline.py ../data/prediction_demo/regression_models/ genomes_retrieved.txt species_taxonomic.txt
```
The final result will be in the file newly_predicted_OGTs.txt, listing each species, the predicted OGT, and the taxonomic model used for the prediction. (Note, these results are not deterministic as the genomes available for each species may change.)

## Model notes
In addition to regression models using all features, models have also been calculated excluding various feature classes (ie. tRNA, rRNA, all RNA, or genome size) which might be absent due to an incomplete organismal genome. These are stored within the subdirectories of data/calculations/prediction/regression_models/ . Therefore, be mindful of your particular particular prediction problem when selecting the regression model directory for the prediction calculation (step 3 above).
