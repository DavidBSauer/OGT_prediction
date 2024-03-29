from sys import argv
import logging
logger = logging.getLogger('feature_calculation')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('feature_calculation.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)
import csv
import multiprocessing as mp
import random
import genomic
import tRNA
import protein
import ORFs
import external_tools
import rRNA
import operator
import os

script, genome_species_file, species_clade_file = argv

#read in the species-clade file. save as dict species->clade
logger.info("trying to open species-taxon file")
infile = open(species_clade_file,'r')
species_clade = {}
for line in infile.readlines():
	working = line.split('\t')
	if not(working[1] == 'None'):
		species_clade[working[0]] = working[1].strip() 
infile.close()    
logger.info("found "+str(len(species_clade.keys()))+" species with taxonomic information")

#read in the genomes file. save a dict of genome->species
logger.info("trying to open genomes file")
infile = open(genome_species_file,'r')
genomes = {}
for line in infile.readlines():
	species = line.split('\t')[1].strip()
	genome = line.split('\t')[0].strip()
	if species in genomes.keys():
		genomes[species].append(genome)
	else:
		genomes[species]=[genome]
infile.close()    
logger.info("found "+str(sum([len(y) for y in genomes.values()]))+" genomes")

#only analyze those with taxonomic description, shuffle
genomes = {species:genomes[species] for species in genomes.keys() if species in species_clade.keys()}
to_analyze = [(genome,species) for species in genomes.keys() for genome in genomes[species]]
random.shuffle(to_analyze)

#make folders for output
if not(os.path.isdir('./output')):	
	os.mkdir('./output')
if not(os.path.isdir('./output/genomes')):	
	os.mkdir('./output/genomes')

for species in [x[1] for x in to_analyze]:
	#make a folder based on the genome name 
	if not(os.path.isdir('./output/genomes/'+species)):	
		os.mkdir('./output/genomes/'+species)

logger.info('Analyzing Genomic features')
#Analyze features per genome
def features_per_genome(inputs):
	(genome_file,species) = inputs
	result = {}
	result['species']=species
	
	#create folders and decompress
	genome_file = external_tools.setup((genome_file,species))	
	
	#calculate genomic features
	logger.info('Working on '+genome_file)
	result['genomic'] = genomic.analysis((genome_file,species))
	
	#calculate tRNA features
	logger.info('Finding tRNA for '+genome_file)
	(tRNA_test,tRNA_seqs) = external_tools.tRNA((genome_file,species))
	if tRNA_test:
		#if tRNAs were predicted, calculate features
		result['tRNA']=tRNA.analysis(tRNA_seqs)

	#calculate rRNA features
	logger.info('Finding rRNAs for '+genome_file)
	#run barrnap using both archaea and bacteria hmm models
	domain_results = external_tools.rRNA((genome_file,species))
	#attempt to predict domain using 16S rRNA sequences
	(pred_test,domain_pred) = external_tools.classify((genome_file,species))
	#calculate rRNA features based on predicted domain
	if pred_test:
		result['rRNA_domain']=domain_pred	
		(rRNA_pred_test,rRNA_pred_seq) = external_tools.rRNA_seq((genome_file,species,domain_pred,'pred'))
		if rRNA_pred_test:		
			result['rRNA_pred']=rRNA.analysis(rRNA_pred_seq)
	#using previous assigned domain
	global species_clade
	domain_assigned = species_clade[species]
	if domain_results[domain_assigned]:
		(rRNA_assigned_test,rRNA_assigned_seq) = external_tools.rRNA_seq((genome_file,species,domain_assigned,'assigned'))
		if rRNA_assigned_test:		
			result['rRNA_assigned']=rRNA.analysis(rRNA_assigned_seq)

	#calculate ORF and proteome features
	logger.info('Identifying and analyzing ORFs for '+genome_file)
	(ORF_test,ORF_seqs) = external_tools.ORF((genome_file,species))
	if ORF_test:		
		t_size = result['genomic']['Total Size']
		result['ORF']=ORFs.analysis(ORF_seqs,t_size)
		(protein_test,protein_seqs) = external_tools.proteins((genome_file,species))
		if protein_test:
			result['protein']=protein.analysis(protein_seqs)
		
	external_tools.cleanup((genome_file,species))
	
	return (species+'/'+genome_file,result)
	
#calculate features in parallel

p = mp.Pool()
results = p.map(features_per_genome, to_analyze)
p.close()
p.join()
'''
#single thread for trouble shooting
results = list(map(features_per_genome,to_analyze))
'''
results = {x[0]:x[1] for x in results}

#write out the results by feature class
for feature_class in ['genomic','tRNA','rRNA_pred','rRNA_assigned','ORF','protein']:
	all_params = []
	for genome in results.keys():
		if feature_class in results[genome].keys():
			all_params = results[genome][feature_class].keys()
			break

	g= open('Genome_'+feature_class+'_features.txt','w')
	g.write('Genome\tspecies\t'+'\t'.join([feature_class+' '+x for x in all_params])+'\n')
	for genome in results.keys():
		if feature_class in results[genome].keys():
			g.write(genome+'\t'+results[genome]['species']+'\t')
			for param in all_params:
				g.write(str(results[genome][feature_class][param])+'\t')
			g.write('\n')
	g.close()

#write out barrnap assignment of domain
h= open('Genome_barrnap_assignment.txt','w')
for genome in results.keys():
	if 'rRNA_domain' in results[genome].keys():
		h.write(genome+'\t'+results[genome]['rRNA_domain']+'\n')
h.close()

logger.info('Exiting normally')
