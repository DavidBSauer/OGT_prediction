# OGT_prediction
Scripts for calculating features from genomic sequences, multiple linear regression of those features to the originating species' Optimal Growth Temperatures (OGT), and prediction of a species' OGT using those linear regression models. 

See: Sauer & Wang. Prediction of Optimal Growth Temperature using only Genome Derived Features (2018) https://doi.org/10.1101/273094

## Requirements
Requires: Python3, BioPython, Numpy, SciPy, Sklearn, MatPlotLib, MatPlotLib-Venn, BCBio-GFF, tqdm, tRNAscan-SE, Bedtools, Barrnap, Prodigal.

In practice, there is a bash script "Ubuntu_setup.bash" that will install everything.

Run as: 
```
./Ubuntu_setup.bash
```

## Usage
### Feature Calculation
Calculate features for each genome
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
(optional) 0. Only when predicting a genome's taxonomy using barrnap identified 16S rRNA sequences, generate a species_taxonomic.txt.
```
python3 genome_species_assignment.py genomes_retrieved.txt Genome_barrnap_assignment.txt
```
1. Calculate correlation of features to OGT and multiple linear regression models.
```
python3 regression_pipeline.py species_trait_file training_testing_file species_taxonomic_file Genome-feature_files
```
training_testing_file is a tab separated file of species annotated as either "train" or "test". For the first run, using NONE will automatically assign 20% of the species to a test set and generate this file for future use.

### Prediction
Predict the OGT of a species (or many species) using the generated multiple linear regression models.
1. Download genomes for species NOT IN a provided list.
```
python3 genome_retriever.py list_of_species_file NOT_IN
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
note: Models for the same taxon in the regression model directory should be cleaned-up so all taxa models are non-redundant prior to use.
