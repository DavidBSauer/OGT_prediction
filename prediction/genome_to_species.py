import logging

def species(genome_features):
	#create lists of features by species
	species_features_list ={}
	all_features =[]
	for genome in genome_features.keys():
		species = genome_features[genome]['species']
		if not(species in species_features_list.keys()):
			species_features_list[species]={}
		curr_features =	[feature for feature in genome_features[genome].keys() if not(feature == 'species')]
		all_features = all_features+curr_features
		for curr_feature in curr_features:	
			if curr_feature in species_features_list[species].keys():#have already seen this feature append
				species_features_list[species][curr_feature].append(genome_features[genome][curr_feature])
			else: #have not seen this feature, create new list
				species_features_list[species][curr_feature]=[genome_features[genome][curr_feature]]

	#in case a feature file is given twice or in pieces
	all_features = list(set(all_features))
	
	#average the lists of feature values for each species
	species_features_average = {}
	for species in species_features_list.keys():
		species_features_average[species] = {}
		for feature in species_features_list[species].keys():
			species_features_average[species][feature] = sum(species_features_list[species][feature])/len(species_features_list[species][feature])

	#record a file of average feature values for all features for each species
	g = open('./species_average_features.txt','w')
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

