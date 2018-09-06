# OGT_prediction
Scripts for calculating features, regression, and prediction of prokaryotic Optimal Growth Temperatures (OGTs)

## Requirements
Requires Python3, BioPython, Numpy, SciPy, Sklearn, MatPlotLib, MatPlotLib-Venn, BCBio-GFF, tqdm, tRNAscan-SE, Bedtools, Barrnap, GeneMarkS

In practice, there is a bash script "Ubuntu_setup.bash" that will install everything except GeneMarkS automatically. 
Run as: 
```
./Ubuntu_setup.bash
```
After installing GeneMarkS, edit the external_tools.txt file with absolute path to gmsn.pl. Then copy external_tools.txt into the feature_calculation and prediction directories.

## Usage
### Feature Calculation
Calculate features for each genome
1. Download genomes for the species of interest.
```
python3 genome_retriever.py list_of_species_file
```
Species file needs to have one per line, with the species name the first thing on the line, where the species name has the form "genus_species" (all lower case, separated by an underscore)
2. Download taxonomic classification for species of interest.
```
python3 clade_retriever.py list_of_species_file
```
3. Calculate feature for each species of interest.
```
python3 feature_calculation_pipeline.py genomes_retrieved.txt species_taxonomic.txt
```
This will produce a series of files with the features for each genome.

### Regression
1. Calculate correlation of features to OGT and multiple linear regression models.
```
python3 regression_pipeline.py species_trait_file training_testing_file species_taxonomic_file Genome-feature_files
```

### Prediction
Predict the OGT of a species (or many species) using the generated multiple linear regression models.
1. Download genomes for the species of interest.
```
python3 genome_retriever.py list_of_species_file
```
Species file needs to have one per line, with the species name the first thing on the line, where the species name has the form "genus_species" (all lower case, separated by an underscore)
2. Download taxonomic classification for species of interest.
```
python3 clade_retriever.py list_of_species_file
```
3. Run the prediction script to predict the OGT of each species.
```
python3 prediction_pipeline.py regression_model_directory/ genomes_retrieved.txt species_taxonomic.txt
```
note: Models for the same taxon in the regression model directory should be cleaned up prior to use.
