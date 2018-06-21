import logging
logging.basicConfig(filename="OGT_calculation.log",filemode="w", level=logging.INFO)

from sys import argv
from sys import exit
import os
import multiprocessing as mp
import genome_analysis
import genome_to_species
from tqdm import tqdm

regression_model_directory = argv[1]
genome_file = argv[2]
species_taxonomy = argv[3]

logging.info('the file of genomic sequences to predict: '+genome_file)
logging.info('the directory of OGT regression models: '+regression_model_directory)
logging.info('file of taxonomy for each genome: '+species_taxonomy)

#load regression models
models = {}
files = [file for file in os.listdir(regression_model_directory) if os.path.isfile(regression_model_directory+file) and file[0] != '.']
for model_file in files:
	logging.info('reading in the '+model_file.split('_')[0]+' linear regression model file: '+regression_model_directory+model_file)
	f = open(regression_model_directory+model_file,'r')
	models[model_file.split('_')[0]] = {line.split('\t')[0].strip():float(line.split('\t')[1].strip()) for line in f.readlines()}
	f.close()		

#build a dictionary of taxonomic information
taxonomic_info={}
infile = open(species_taxonomy,'r')
ranks = [x.strip() for x in infile.readline().split('\t')[1:]]
for line in infile.readlines():
	species = line.split('\t')[0]
	if not(species in taxonomic_info.keys()):
		taxonomic_info[species] ={}
	taxons = [x.strip() for x in line.split('\t')[1:]]
	for x in range(0,len(ranks),1):
		if not(taxons[x] =='None'):
			taxonomic_info[species][ranks[x]]=taxons[x]
infile.close()

#calculate genomic features for each genome
to_analyze ={}
f =open(genome_file,'r')
for line in f.readlines():
	to_analyze[line.split('\t')[0].strip()]=line.split('\t')[1].strip()
f.close()
logging.info('initial number of species to be predicted: '+str(len(list(set(to_analyze.values())))))

logging.info('initial number of genomes to be analyzed: '+str(len(to_analyze.keys())))
to_analyze_tax = {species:taxonomic_info[species]['superkingdom'] for species in taxonomic_info.keys() if (('superkingdom' in taxonomic_info[species].keys()) and (taxonomic_info[species]['superkingdom'] in ['Bacteria','Archaea']) and (species in to_analyze.values()))} #analyze only bacteria and archaea, and require that the superkingdom be assigned
logging.info('number of species to be predicted: '+str(len(list(set(to_analyze_tax.keys())))))
genome_analysis.setup(to_analyze_tax)
logging.info('number of genomes to be analyzed: '+str(len([genome for genome in to_analyze.keys() if to_analyze[genome] in to_analyze_tax.keys()])))
genome_features = genome_analysis.many_genomes(to_analyze)

#reduce genome features to species features
species_features = genome_to_species.species(genome_features)
logging.info('number of species with features calculated: '+str(len(species_features)))

#calculate OGT for each species from previous feature linear regressions
newly_predicted_OGTs = {}
print 'calculating OGTs'
for species in tqdm(species_features.keys(),unit='species'):
	for rank in ['superkingdom','phylum','class']: #will continually overwrite the value with a progressively more taxon specific prediction
		if rank in taxonomic_info[species].keys():
			if taxonomic_info[species][rank] in models.keys():
				regress_features = models[taxonomic_info[species][rank]].keys()
				regress_features_less = [feature for feature in regress_features if not(feature == 'intercept')]
				if set(regress_features_less).issubset(set(species_features[species].keys())):
					newly_predicted_OGTs[species]={}
					pred_ogt = []
					for feature in regress_features_less:
						pred_ogt.append(models[taxonomic_info[species][rank]][feature]*species_features[species][feature])
					pred_ogt.append(models[taxonomic_info[species][rank]]['intercept'])
					pred_ogt = sum(pred_ogt)
					newly_predicted_OGTs[species]['OGT'] = pred_ogt
					newly_predicted_OGTs[species]['model'] = taxonomic_info[species][rank]

g = open('newly_predicted_OGTs.txt','w')
for species in newly_predicted_OGTs.keys():
	g.write(species+'\t'+str(newly_predicted_OGTs[species]['OGT'])+'\t'+newly_predicted_OGT[species]['model']+'\n')
g.close()
logging.info('number of species with a newly predicted OGT: '+str(len(newly_predicted_OGTs.keys())))

