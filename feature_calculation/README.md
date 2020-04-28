# Feature Calculation
Calculates features for a provided list of species.

## Usage
### Feature Calculation
Calculate features for each genome. Start within the feature_calculation directory.

1. Download genomes for species IN a provided list.
```
python3 genome_retriever.py list_of_species_file IN
```
list_of_species_file needs to have one species per line, with the species name the first thing on the line, where the species name has the form "genus_species" (all lower case, separated by an underscore). The release is the release number of Ensembl Bacteria to use. Release 40 was used for Sauer & Wang (2019).

2. Download taxonomic classification for species of interest.
```
python3 clade_retriever.py species_retrieved.txt your_email_address@awesome.com
```

3. Calculate feature for each species of interest.
```
python3 feature_calculation_pipeline.py genomes_retrieved.txt species_taxonomic.txt
```
This will produce a series of files with the features for each genome. The feature files will be titled "Genome_XXX_features.txt" where XXX is the feature class (genomic, tRNA, rRNA, ORF, protein).
