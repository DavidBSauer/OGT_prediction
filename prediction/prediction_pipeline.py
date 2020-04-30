import logging
logger = logging.getLogger('prediction')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('prediction.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)

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

logger.info('the directory of OGT regression models: '+regression_model_directory)
logger.info('the file of genomic sequences to predict: '+genome_file)
logger.info('file of taxonomy for each genome: '+species_taxonomy)
r2_limit = 0.5
logger.info('only using regression models with an R2 >= '+str(r2_limit))

#load regression models
models = {}
files = [file for file in os.listdir(regression_model_directory) if os.path.isfile(regression_model_directory+file) and file[0] != '.']
for model_file in files:
	rank = model_file.split('-')[0]
	if not(rank in models.keys()):
		models[rank]={}
	clade = model_file.split('-')[1]
	models[rank][clade] = {}
	logger.info('reading in the '+rank+'-'+clade+' linear regression model file: '+regression_model_directory+model_file)
	f = open(regression_model_directory+model_file,'r')
	lines = f.readlines()
	f.close()
	models[rank][clade]['R2']=float(lines[0].split('|')[0].split('=')[1])
	models[rank][clade]['model'] = {line.split('\t')[0].strip():float(line.split('\t')[1].strip()) for line in lines[1:]}
	

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

#read in the genomes available to be analyzed
to_analyze ={}
f =open(genome_file,'r')
count = 0
for line in f.readlines():
	#to_analyze[line.split('\t')[0].strip()]=line.split('\t')[1].strip()
	species = line.split('\t')[1].strip()
	genome = line.split('\t')[0].strip()
	if species in to_analyze.keys():
		to_analyze[species].append(genome)
	else:
		to_analyze[species]=[genome]
	count = count+1
f.close()
logger.info('number of lines in the genome-species file: '+str(count))

#find the genomes and species to be analyzed
logger.info('initial number of species to be predicted: '+str(len(to_analyze.keys())))
logger.info('initial number of genomes to be analyzed: '+str(sum([len(y) for y in to_analyze.values()])))
useful_tax_info = {species:taxonomic_info[species] for species in taxonomic_info.keys() if 'superkingdom' in taxonomic_info[species].keys()} #keep only species with an assigned superkingdom
useful_tax_info = {species:useful_tax_info[species] for species in useful_tax_info.keys() if useful_tax_info[species]['superkingdom'] in ['Bacteria','Archaea']} #keep only species which are archaea or bacteria
to_analyze = {species:to_analyze[species] for species in to_analyze.keys() if species in useful_tax_info.keys()} #analyze only those genomes with available species taxonomic information
logger.info('number of species to be predicted: '+str(len(to_analyze.keys())))
logger.info('number of genomes to be analyzed: '+str(sum([len(y) for y in to_analyze.values()])))

#calculate genome features
genome_analysis.setup(useful_tax_info)
genome_features = genome_analysis.many_genomes(to_analyze)
#reduce genome features to species features
species_features = genome_to_species.species(genome_features)
logger.info('number of species with features calculated: '+str(len(species_features.keys())))

#calculate OGT for each species from previous feature linear regressions
newly_predicted_OGTs = {}
print('calculating OGTs')
for species in tqdm(species_features.keys(),unit='species'):
	for rank in ['domain','superkingdom','phylum','class','order','family']: #will continually overwrite the value with a progressively more taxon specific prediction
		if ((rank in list(taxonomic_info[species].keys())) and (rank in models.keys())):
			if taxonomic_info[species][rank] in list(models[rank].keys()):
				if models[rank][taxonomic_info[species][rank]]['R2'] >= r2_limit:
					#logger.info('Testing '+species+' against the model '+rank+'-'+taxonomic_info[species][rank])
					regress_features = models[rank][taxonomic_info[species][rank]]['model'].keys()
					regress_features_less = [feature for feature in regress_features if not(feature == 'intercept')]
					if set(regress_features_less).issubset(set(species_features[species].keys())):
						logger.info(species+' has all the features for the model '+rank+'-'+taxonomic_info[species][rank]+'. predicting OGT')
						newly_predicted_OGTs[species]={}
						pred_ogt = []
						for feature in regress_features_less:
							pred_ogt.append(models[rank][taxonomic_info[species][rank]]['model'][feature]*species_features[species][feature])
						pred_ogt.append(models[rank][taxonomic_info[species][rank]]['model']['intercept'])
						pred_ogt = sum(pred_ogt)
						newly_predicted_OGTs[species]['OGT'] = pred_ogt
						newly_predicted_OGTs[species]['model'] = rank+'-'+taxonomic_info[species][rank]
					#else:
						#logger.info(species+' does not have all the features to use the model '+rank+'-'+taxonomic_info[species][rank])
				#else:
					#logger.info('While the model '+rank+'-'+taxonomic_info[species][rank]+' is specific to '+species+', the R2 of the model is < 0.5. Therefore skipping this model.')



g = open('newly_predicted_OGTs.txt','w')
for species in newly_predicted_OGTs.keys():
	g.write(species+'\t'+str(newly_predicted_OGTs[species]['OGT'])+'\t'+newly_predicted_OGTs[species]['model']+'\n')
g.close()
logger.info('number of species with a newly predicted OGT: '+str(len(newly_predicted_OGTs.keys())))

