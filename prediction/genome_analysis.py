from sys import argv
import logging
logger = logging.getLogger('prediction')
import csv
import multiprocessing as mp
import genomic
import tRNA
import protein
import ORFs
import external_tools
import rRNA
import operator
import random
import os
from tqdm import tqdm

species_clade = {}
def setup(sc_input):
	global species_clade
	species_clade = sc_input

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
	domain_results = external_tools.rRNA((genome_file,species))
	#using previous assigned domain
	global species_clade
	domain_assigned = species_clade[species]['superkingdom']
	if domain_results[domain_assigned]:
		(rRNA_assigned_test,rRNA_assigned_seq) = external_tools.rRNA_seq((genome_file,species,domain_assigned,'assigned'))
		if rRNA_assigned_test:		
			result['rRNA']=rRNA.analysis(rRNA_assigned_seq)

	#calculate ORF and proteome features
	logger.info('Identifying and analyzing ORFs for '+genome_file)
	(ORF_test,ORF_seqs) = external_tools.prodigal((genome_file,species))
	if ORF_test:		
		t_size = result['genomic']['Total Size']
		result['ORF']= ORFs.analysis(ORF_seqs,t_size)
		(protein_test,protein_seqs) = external_tools.proteins((genome_file,species))
		if protein_test:
			result['protein']=protein.analysis(protein_seqs)
		
	external_tools.cleanup((genome_file,species))
	
	return (genome_file,result)

def many_genomes(genomes):
	to_analyze = [(genome,genomes[genome]) for genome in genomes]
	random.shuffle(to_analyze)
	print('analyzing genomes')
	#make folders for output
	if not(os.path.isdir('./output')):	
		os.mkdir('./output')
	if not(os.path.isdir('./output/genomes')):	
		os.mkdir('./output/genomes')

	for species in list(set([x[1] for x in to_analyze])):
		#make a folder based on the genome name 
		if not(os.path.isdir('./output/genomes/'+species)):	
			os.mkdir('./output/genomes/'+species)
	
	#calculate features in parallel
	p = mp.Pool()
	results = p.map(features_per_genome, to_analyze)
	p.close()
	p.join()
	'''
	#calculate single thread for troubleshooting
	results =[]
	for x in tqdm(to_analyze,unit='genome'):
		results.append(features_per_genome(x))
	'''
	results = {x[0]:x[1] for x in results}
	
	remod_results ={genome:{} for genome in results.keys()}
	#write out the results by feature class
	for feature_class in ['genomic','tRNA','rRNA','ORF','protein']:
		all_params = []
		for genome in results.keys():
			if feature_class in results[genome].keys():
				all_params = results[genome][feature_class].keys()
				break

		g= open('./output/Genome_'+feature_class+'_features.txt','w')
		g.write('Genome\tspecies\t'+'\t'.join([feature_class+' '+x for x in all_params])+'\n')
		for genome in results.keys():
			if feature_class in results[genome].keys():
				g.write(genome+'\t'+results[genome]['species']+'\t')
				for param in all_params:
					g.write(str(results[genome][feature_class][param])+'\t')
					remod_results[genome][feature_class+' '+param] = results[genome][feature_class][param]
				g.write('\n')
		g.close()
	for genome in results.keys():
		remod_results[genome]['species'] = results[genome]['species']
	return remod_results
