import logging

#calculate histograms of OGTS vs genomes and species

def calc(genome_species,species_trait,species_clade,unit):
	#calculate genome_redundancy vs species count
	species_count = {}
	genomes = genome_species.keys()
	for genome in genomes:
		if genome_species[genome] in species_count.keys():
			species_count[genome_species[genome]] = species_count[genome_species[genome]]+1
		else:
			species_count[genome_species[genome]] = 1
	genome_redundancy = {x:0 for x in species_count.values()}
	for x in species_count.keys():
		genome_redundancy[species_count[x]] = genome_redundancy[species_count[x]]+1
	g= open('./files/genome_redundancy.txt','w')
	g.write('Genome_redundancy\tSpecies_count\n')
	keys = genome_redundancy.keys()
	keys.sort()
	for x in keys:
		g.write(str(x)+'\t'+str(genome_redundancy[x])+'\n')
	g.close()

	#write a file of the species-number of genomes
	g = open('./files/species_genomes.txt','w')
	g.write('Species\tGenomes\n')
	for species in species_count.keys():
		g.write(species+'\t'+str(species_count[species])+'\n')
	g.close()

	#create a list of all species
	all_species = list(set([genome_species[genome] for genome in genomes]))
	all_species = [species for species in all_species if species in species_trait.keys()]
	
	#automatically find bounds of range
	min_bound = int(0.95*min(species_trait.values()))
	max_bound = int(1.05*max(species_trait.values()))
	
	#create histogram with 25 steps
	all_buckets = {q:0 for q in range(min_bound,max_bound,int(max_bound-min_bound)/25)}
	for clade in ['all','Bacteria','Archaea']:
		g = open('./files/histogram_'+clade+'.txt','w')
		buckets = {q:0 for q in range(min_bound,max_bound,int(max_bound-min_bound)/25)}
		for species in all_species:
			if species in species_clade['superkingdom'][clade]:
				for z in range(min_bound,max_bound,int(max_bound-min_bound)/25):
					if z+5 > species_trait[species] >= z :
						buckets[z]=buckets[z]+1
						all_buckets[z]=all_buckets[z]+1
		for q in range(min_bound,max_bound,int(max_bound-min_bound)/25):
			g.write(str(q)+' to '+str(q+int(max_bound-min_bound)/25)+' count: '+str(buckets[q])+'\n')
		g.close()

	#log the minimum and maximum OGTs
	minimum = min([species_trait[species] for species in all_species])
	maximum = max([species_trait[species] for species in all_species])
	logging.info('The mimimum reported '+unit+' in the dataset: '+str(minimum))
	logging.info('The maximum reported '+unit+' in the dataset: '+str(maximum))
	
	#record the number in each superkingdom
	bacteria = len([species for species in all_species if species in species_clade['superkingdom']['Bacteria']])
	archaea = len([species for species in all_species if species in species_clade['superkingdom']['Archaea']])
	
	logging.info('The number of species in the dataset: '+str(len(all_species)))
	logging.info('The number of bacteria in the dataset: '+str(bacteria))
	logging.info('The number of archaea in the dataset: '+str(archaea))
