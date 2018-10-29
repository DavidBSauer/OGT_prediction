# Regression
Calculates feature correlaton to OGT and multiple linear regression models for species' features to OGT.

## Usage
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
