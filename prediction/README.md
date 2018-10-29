# Prediction 
Predict the OGT for a provided set of species, given a set of multiple linear regression models.

## Usage
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
