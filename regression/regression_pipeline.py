import logging
logging.basicConfig(filename='regression_pipeline.log', level=logging.INFO)
from sys import argv
from sys import version
import genome_to_species as feature_average
import feature_regression_parallel as feature_regression
import multi_lin_regression
import random
import feature_assignment_venn as venn
import csv
import histogram_calc as histogram
import cross_corr

logging.info("python version: "+version)
#run as python regression_pipeline species_OGT_file train_test_file species_taxon_file genome_features tRNA_features rRNA_features ORF_features Protein_features

#whole pipeline
species_OGT_file = argv[1]
logging.info('Species-OGT file: '+species_OGT_file)
train_test_species_file = argv[2] #may be NONE if not using previous generated train-test list.
logging.info('Training-test file: '+train_test_species_file)
species_taxon_file = argv[3]
logging.info('Species-taxon file: '+species_taxon_file)
feature_files = argv[4:]
for file in feature_files:
	logging.info('found species-feature file of: '+file)

#parse species-OGT file. creates dict of species -> OGT
infile = open(species_OGT_file,'r')
reader = csv.reader(infile,delimiter='\t')
species_OGT = dict((str(rows[0]),float(rows[1])) for rows in reader)
infile.close()  
logging.info("found "+str(len(species_OGT.keys()))+" species-OGT pairs")
unit = 'OGT'

#given a file of taxonomic assignment for each species create a dictionary of the same form as previously: taxon level -> clade -> list of species for that clade
all_species = []
infile = open(species_taxon_file,'r')
levels = [x.strip() for x in infile.readline().split('\t')[1:]]
species_rank_dicts = {level:{} for level in levels}
for line in infile.readlines():
	species = line.split('\t')[0]
	if species in species_OGT.keys():
		all_species.append(species)
		values = [x.strip() for x in line.split('\t')[1:]]
		for x in range(0,len(levels),1):
			if not(values[x] == 'None'):
				if values[x] in species_rank_dicts[levels[x]].keys():
					species_rank_dicts[levels[x]][values[x]].append(species)
				else:
					species_rank_dicts[levels[x]][values[x]] = [species]
infile.close()
for level in levels:
	logging.info("for rank of "+level+" found "+str(len(species_rank_dicts[level].keys()))+" taxonomic groups")
species_rank_dicts['superkingdom']['all']=list(set(all_species))
#remove empty taxons
species_rank_dicts['phylum'].pop('', None)
species_rank_dicts['class'].pop('', None)
logging.info("found "+str(len(species_rank_dicts['superkingdom']['all']))+" total species with taxonomic assignment")
	
#parse/create train-test file
testing = []
#if a train-test list does not exist, generate one. otherwise use the previously created list.
if train_test_species_file == 'NONE':
	logging.info('Assigning new Train-test list. Saving as train_test_list.txt')
	f = open('train_test_list.txt','w')
	for species in species_OGT.keys():
		if random.random() >0.8:
			testing.append(species)
			f.write(species+'\ttest\n')
		else:
			f.write(species+'\ttrain\n')
	f.close()
else:
	infile = open(train_test_species_file,'r')
	testing = [line.split('\t')[0].strip() for line in infile.readlines() if line.split('\t')[1].strip() == 'test']
	infile.close()

#average redundant genome-features -> return species->feature->value dictionary. calculate genome and species histograms
logging.info('Calculate species average features from all genomes')
species_features,features,genome_species = feature_average.species(feature_files)
logging.info('Number of species: '+str(len(species_features.keys())))

#remove those species without a OGT recorded
species_features = {species:species_features[species] for species in species_features.keys() if species in species_OGT.keys()}

#calculate histogram of OGTs
logging.info('Calculate histograms of OGTs')
histogram.calc(genome_species,species_OGT,species_rank_dicts,unit)

#calculate venn diagram of data
logging.info('Calculate data venn diagram')
venn.calc(species_features,species_OGT,species_rank_dicts)

analysis_sets = [['genomic'],['genomic','tRNA'],['genomic','rRNA'],['genomic','tRNA','rRNA'],['genomic','tRNA','rRNA','ORF'],['genomic','tRNA','rRNA','ORF','protein']]
#calculate regression for all superkingdoms

logging.info('Calculate multiple linear regression of all data')
for clade in ['all']:
	rank ='superkingdom'
	#calculate rs
	logging.info('Calculate feature r-values to OGT for domain: '+clade)
	rvalues = feature_regression.rs(features,species_features,species_OGT,rank,clade,species_rank_dicts[rank][clade],unit)

	#calculate feature cross-correlation R values with figure
	logging.info('Calculate feature cross correlation |r| values')
	cross_corr.calc(rvalues,features,species_features,analysis_sets[-1],clade)

	#for all species, calculated with increasing feature classes
	for analysis in analysis_sets:
		logging.info('Running multiple linear regression for '+clade+' species and feature sets of '+'+'.join(analysis))
		multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis),species_rank_dicts[rank][clade],analysis,features,species_features,species_OGT,testing,rvalues,unit)
	
	#limit OGT range to species with OGT >= 25C
	analysis = analysis_sets[-1]
	logging.info('Running multiple linear regression for '+rank+'-'+clade+' species with OGT >25C')
	valid_species = [species for species in species_rank_dicts[rank][clade] if species_OGT[species]>25]	
	multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis)+'_>25',valid_species,analysis,features,species_features,species_OGT,testing,rvalues,unit)

	#exclude genome size
	logging.info('Running multiple linear regression for '+rank+'-'+clade+' species excluding genome size')
	del rvalues['genomic Total Size']
	multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis)+'_ex_genome_size',species_rank_dicts[rank][clade],analysis,features,species_features,species_OGT,testing,rvalues,unit)

logging.info('Calculate multiple linear regression of superkingdom specific data')
for clade in ['Archaea','Bacteria']:
	rank ='superkingdom'
	#calculate rs
	logging.info('Calculate feature r-values to OGT for domain: '+clade)
	rvalues = feature_regression.rs(features,species_features,species_OGT,rank,clade,species_rank_dicts[rank][clade],unit)

	#for specific superkingdoms, calculate regression using all features classes
	analysis = analysis_sets[-1]
	logging.info('Running multiple linear regression for '+rank+'-'+clade+' species and feature sets of '+'+'.join(analysis))
	multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis),species_rank_dicts[rank][clade],analysis,features,species_features,species_OGT,testing,rvalues,unit)

	#limit OGT range to species with OGT >= 25C
	analysis = analysis_sets[-1]
	logging.info('Running multiple linear regression for '+rank+'-'+clade+' species with OGT >25C')
	valid_species = [species for species in species_rank_dicts[rank][clade] if species_OGT[species]>25]	
	multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis)+'_>25',valid_species,analysis,features,species_features,species_OGT,testing,rvalues,unit)

	#exclude genome size
	logging.info('Running multiple linear regression for '+clade+' species excluding genome size')
	del rvalues['genomic Total Size']
	multi_lin_regression.regress(rank+'_'+clade+'_'+'+'.join(analysis)+'_ex_genome_size',species_rank_dicts[rank][clade],analysis,features,species_features,species_OGT,testing,rvalues,unit)


#run regression by clade
for rank in [level for level in species_rank_dicts.keys() if not(level == 'superkingdom')]:
	logging.info('Calculate multiple linear regression of '+rank+' specific data')
	for clade in species_rank_dicts[rank].keys():
		valid_species = [species for species in species_rank_dicts[rank][clade] if species in species_features.keys()] #find those species with features calculated	
		if len(valid_species)>=100: #only run regression if >= 50 species
			logging.info('Calculate feature r-values for '+rank+'-'+clade)
			rvalues = feature_regression.rs(features,species_features,species_OGT,rank,clade,species_rank_dicts[rank][clade],unit)
			if len([x for x in rvalues.keys() if abs(rvalues[x])>0.3])>0: #run regression only if there are correlated features
				analysis = analysis_sets[-1] #calculate regression using all features
				logging.info('Running multiple linear regression for '+rank+'-'+clade+' species and feature sets of '+'+'.join(analysis))
				multi_lin_regression.regress(rank+'_'+clade+'_all_features',species_rank_dicts[rank][clade],analysis,features,species_features,species_OGT,testing,rvalues,unit)
		
