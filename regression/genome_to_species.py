import logging
logger = logging.getLogger('regression')

#take in file names for features and genomes_species -> return dictionary of species-feature and list of all features
def species(feature_files):
	all_features=[]
	species_features_list = {}
	genomes_species ={}
	for feature_file in feature_files:
		#read in the feature file
		f= open(feature_file,'r')
		features = [x.strip() for x in f.readline().split('\t')[2:]]
		
		#make rRNA generic
		features = [x.replace('rRNA_assigned','rRNA') for x in features]
		features = [x.replace('rRNA_pred','rRNA') for x in features]
		
		all_features = all_features+features #add features to list of all features
		
		#working through each line of the feature file, record the values
		for line in f.readlines():
			working_line = line.strip() #remove any dangling whitespace
			working_line = [x.strip() for x in working_line.split('\t')]
			data = [float(x) for x in working_line[2:]] #convert values to floats
			genome = working_line[0]
			species = working_line[1]
			genomes_species[genome] = species #map each genome to a species
			#if the species has already been seen previously
			if species in species_features_list.keys():
				for z in range(0,len(features),1):
					if features[z] in species_features_list[species]: #if the feature has already been seen previously
						species_features_list[species][features[z]].append(data[z])
					else: #if the feature has not been seen previously
						species_features_list[species][features[z]]=[data[z]]
			else: #if the species has not been seen previously
				species_features_list[species] = {}
				for z in range(0,len(features),1):
					species_features_list[species][features[z]]=[data[z]]
	
	#in case a feature file is given twice or in pieces
	all_features = list(set(all_features))
	
	#average the lists of feature values for each species
	species_features_average = {}
	for species in species_features_list.keys():
		species_features_average[species] = {}
		for feature in species_features_list[species].keys():
			species_features_average[species][feature] = sum(species_features_list[species][feature])/len(species_features_list[species][feature])

	#record a file of average feature values for all features for each species
	g = open('./files/species_average_features.txt','w')
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
	return (species_features_average,all_features,genomes_species)

