# OGT_prediction
Scripts for calculating features from genomic sequences, multiple linear regression of those features to the originating species' Optimal Growth Temperatures (OGT), and prediction of a species' OGT using those linear regression models. 

See: Sauer & Wang. Prediction of Optimal Growth Temperature using only Genome Derived Features (2018) https://doi.org/10.1101/273094

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
./Ubuntu_setup.bash
```
If you're on another OS, install all the python packages. Then install bedtools, tRNAscan-SE, barrnap, and prodigal; and all their dependencies. Create a file called external_tools.txt listing bedtools, tRNAscan-SE, barrnap, and prodigal tab separated from the absolute path for each executable. Copy this external_tools.txt file into the feature_calculation and prediction directories.

## Usage
### Feature Calculation
Calculate features for each genome. Start within the feature_calculation directory.

1. Download genomes for species IN a provided list.
```
python3 genome_retriever.py list_of_species_file IN
```
list_of_species_file needs to have one species per line, with the species name the first thing on the line, where the species name has the form "genus_species" (all lower case, separated by an underscore).

2. Download taxonomic classification for species of interest.
```
python3 clade_retriever.py species_retrieved.txt your_email_address@awesome.com
```

3. Calculate feature for each species of interest.
```
python3 feature_calculation_pipeline.py genomes_retrieved.txt species_taxonomic.txt
```
This will produce a series of files with the features for each genome. The feature files will be titled "Genome_XXX_features.txt" where XXX is the feature class (genomic, tRNA, rRNA, ORF, protein).

### Regression
Start within the regression directory.

(optional) 0. Only when predicting a genome's taxonomy using barrnap identified 16S rRNA sequences, generate barrnap_species_taxonomic.txt.
```
python3 genome_species_assignment.py genomes_retrieved.txt Genome_barrnap_assignment.txt
```

1. Calculate correlation of features to OGT and multiple linear regression models.
```
python3 regression_pipeline.py species_trait_file training_testing_file species_taxonomic_file Genome-feature_files
```
training_testing_file is a list of species annotated as either "train" or "test", tab separated. For the first run, using NONE will automatically assign 20% of the species to a test set and generate this file for future use.

### Prediction
Predict the OGT of a species (or many species) using the generated multiple linear regression models. Start within the prediction directory.

1. Download genomes for species IN or NOT_IN a provided list.
```
python3 genome_retriever.py list_of_species_file IN/NOT_IN
```
Species file needs to have one per line, with the species name the first thing on the line, where the species name has the form "genus_species" (all lower case, separated by an underscore).

2. Download taxonomic classification for species of interest.
```
python3 clade_retriever.py species_retrieved.txt your_email_address@awesome.com
```

3. Run the prediction script to predict the OGT of each species.
```
python3 prediction_pipeline.py regression_model_directory/ genomes_retrieved.txt species_taxonomic.txt
```
note: Models for the same taxon in the regression model directory should be cleaned-up so all taxa models are non-redundant prior to use. Each model file should be named rank-taxon-xxx. Each model file is a list features-coefficients pairs, tab separated.

The final result will be in the file newly_predicted_OGTs.txt, listing each species, the predicted OGT, and the taxonomic model used for the prediction.

## Prediction Demonstration
This uses the previously computed regression models to predict OGTs for a few species. Start within the prediction directory.

1. Download genomes species IN the provided list.
```
python3 genome_retriever.py ../data/prediction_demo/species.txt IN
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
