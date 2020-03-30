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

## OGT Data
To calculate regressions, a species_trait file (tab delimited) is required. The species_OGT file used for calculating regression in Sauer & Wang (2019) and found here, was originally published in Sauer et al. (2015). The species_OGT file is named "all_merged_12_10_2012.txt" and found within http://www.med.nyu.edu/skirball-lab/dwanglab/files/thermostability_v1.1.tar.gz. If new regression are calcuated using that species-OGT file, please cite: Sauer et al. Rapid Bioinformatic Identification of Thermostabilizing Mutations. Biophysical Journal (2015) https://doi.org/10.1016/j.bpj.2015.07.026.
