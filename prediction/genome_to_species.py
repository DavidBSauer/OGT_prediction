import logging
logger = logging.getLogger('prediction')

def species(genome_features):
	all_features=[]
	species_features_list = {}
	#read the genome->feature dict
	for genome in genome_features.keys():
		features = list(genome_features[genome].keys())
		features = [x for x in features if not(x == 'species')]
		#make rRNA generic
		features = [x.replace('rRNA_assigned','rRNA') for x in features]
		features = [x.replace('rRNA_pred','rRNA') for x in features]
		all_features = all_features+features
		species = genome_features[genome]['species']
		if species in species_features_list.keys():
			for feature in features:
				if feature in species_features_list[species]: #if the feature has already been seen previously
					species_features_list[species][feature].append(genome_features[genome][feature])
				else: #if the feature has not been seen previously
					species_features_list[species][feature]=[genome_features[genome][feature]]
		else: #if the species has not been seen previously
			species_features_list[species] = {}
			for feature in features:
				species_features_list[species][feature]=[genome_features[genome][feature]]

	#in case a feature file is given twice or in pieces
	all_features = list(set(all_features))
	
	#average the lists of feature values for each species
	species_features_average = {}
	for species in species_features_list.keys():
		species_features_average[species] = {}
		for feature in species_features_list[species].keys():
			species_features_average[species][feature] = sum(species_features_list[species][feature])/len(species_features_list[species][feature])

	#record a file of average feature values for all features for each species
	g = open('./output/species_average_features.txt','w')
	g.write(' \t'+'\t'.join(all_features)+'\n')
	for species in species_features_average.keys():
		g.write(species)
		for feature in all_features:
			g.write('\t')
			if feature in species_features_average[species].keys():
				g.write(str(species_features_average[species][feature]))
			else:
				g.write('NONE')
		g.write('\n')
	g.close()
	return species_features_average

