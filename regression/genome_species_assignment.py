#convert a barrnap list of genome-superkingdom into a conventional list of species-superkingdom
import logging
logger = logging.getLogger('taxon_assignment')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('taxon_assignment.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)
from sys import argv

species_genome_file, genome_assignment_file = argv[1:3]
logger.info('Genome-species file: '+species_genome_file)
logger.info('Genome-barrnap assignment file: '+genome_assignment_file)

f = open(species_genome_file,'r')
genome_species = {'.'.join(line.strip().split('\t')[0].split('.')[:-1]):line.strip().split('\t')[1] for line in f.readlines()}
f.close()
logger.info('Number of genomes with an assigned species: '+str(len(genome_species.keys())))

f = open(genome_assignment_file,'r')
genome_assignment = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in f.readlines()}
f.close()
logger.info('Number of genomes with a barrnap assigned taxon: '+str(len(genome_assignment.keys())))

common_genomes = [genome for genome in genome_assignment.keys() if genome in genome_species.keys()]
logger.info('Number of genomes in common on both lists: '+str(len(common_genomes)))

species_taxon_lists = {}
for genome in common_genomes:
	if genome_species[genome] in species_taxon_lists.keys():
		species_taxon_lists[genome_species[genome]].append(genome_assignment[genome])
	else:
		species_taxon_lists[genome_species[genome]]=[genome_assignment[genome]]

species_taxon = {}
for species in species_taxon_lists.keys():
	assignment = 'Bacteria' #default to bacteria
	if species_taxon_lists[species].count('Archaea')>species_taxon_lists[species].count('Bacteria'):
		assignment = 'Archaea'
	species_taxon[species] = assignment
logger.info('Number of species with an assigned taxon: '+str(len(species_taxon.keys())))

f = open('barrnap_species_taxonomic.txt','w')
f.write('species\tsuperkingdom\n')
for species in species_taxon.keys():
	f.write(species+'\t'+species_taxon[species]+'\n')
f.close()
